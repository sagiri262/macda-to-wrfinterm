"""MACDA v2.0 conversion driver."""

# 本模块保留 MACDA 转 WPS intermediate 的主流程，并提供 MCD 共用的 CSV 变量表读取器。
# MACDA 路径依次执行时间核验、坐标方向统一、单位/物理转换、sigma 到等压面插值、水平插值。
# MCD 的混合压力不能直接套用 MACDA sigma 算法，应由 mcd_pipeline 单独准备和导出。
# 文件写出后的回读检查只验证中间文件契约，不等于 real.exe 已接收或生成了正确初始场。

# 推迟求值类型注解，允许在签名里使用 list[Path] 和联合类型等标注。
from __future__ import annotations

# CSV 描述字段映射，dataclass 保存一行映射，Path/Iterable 用于路径和时间索引接口。
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# 配置对象提供目标网格、压力层、时间策略和可选字段开关；底层模块负责各自的具体算法。
from macda2wrf.config import MacdaConfig
from macda2wrf.dust import (columnOpticalDepth, combineDustTracers,
                            massExtinctionCoefficient, opticalDepthScaling,
                            tracerFromDustMass, validateDustPartition)
from macda2wrf.grid import RegularLatLonGrid, horizontal_interp_2d, orient_lat_lon
from macda2wrf.macda_reader import MacdaReader
from macda2wrf.mars_time import (ANCHORED_HDATE_STRATEGIES, LabelledRecord,
                                 MacdaRecord, MARSWRF_FIXED_YEAR_SOLS,
                                 build_label_sequence, short_mars_years_left, validate_macda_time)
from macda2wrf.vertical import interp_sigma_to_pressure
from macda2wrf.wrf_intermediate import (IntermediateRecord, WrfIntermediateWriter,
                                        XLVL_SURFACE, validate_intermediate_file)

DUST_BIN1_TRANSFORM = "dust_to_bin1_tracer"
DUST_BIN2_TRANSFORM = "dust_to_bin2_tracer"
DUST_TRANSFORMS = {DUST_BIN1_TRANSFORM, DUST_BIN2_TRANSFORM}


@dataclass(frozen=True)
class FieldSpec:
    # frozen=True 禁止加载后修改映射对象，防止运行途中各时刻采用不同的字段定义。
    # src_v 为源变量名或约定占位名，aim_v 为 WPS 字段名，units 写入中间文件头。
    src_v: str
    aim_v: str
    units: str
    # kind 决定二维/三维及垂直坐标处理分支；transform 决定字段的数值变换。
    kind: str
    transform: str
    # required 表示必需字段，enabled 控制该行是否启用；desc 为 WPS 字段说明。
    required: bool
    enabled: bool
    desc: str


class MacdaConverter:
    def __init__(self, cfg: MacdaConfig):
        # 只保存配置并建立目标网格描述，此时尚未读取整份 MACDA 数据或写文件。
        self.cfg = cfg
        # 起止经纬度与格点数共同决定等间距网格；火星半径以 km 传入 WPS 网格头。
        self.grid = RegularLatLonGrid(lat_start=cfg.lat_start,
            lat_end=cfg.lat_end, nlat=cfg.nlat,
            lon_start=cfg.lon_start, lon_end=cfg.lon_end,
            nlon=cfg.nlon, radius_km=cfg.mars_radius_km,)
        # 全部字段先读入；enabled/required 与 emit_optional_fields 在实际输出时共同筛选。
        self.fields = read_field_table(cfg.variable_table)
        validateDustFieldTable(self.fields)

    def run(self, time_indices: Iterable[int] | None = None,
            max_times_override: int | None = None,
            dry_run: bool = False) -> list[Path]:
        # 返回成功写出的路径列表。输入可显式指定时间索引；dry_run 只核对和打印处理计划。
        written: list[Path] = []
        with MacdaReader(self.cfg.input_file) as reader:
            # 先读小体积坐标；lev 是 MACDA sigma 层，三维压力在每个时刻由 sigma*ps 得到。
            lat = reader.read("lat")
            lon = reader.read("lon")
            sigma = reader.read("lev")
            time_values = reader.read("time")
            ntimes = reader.size("time")
            anchor_index = self._anchor_index(ntimes)
            # 显式索引优先，否则从配置选闭区间；命令行 max 覆盖存在时不再先截断配置上限。
            indices = self._select_indices(ntimes, anchor_index,
                time_indices, max_times_override)
            labelled = self._build_labels(reader, indices, time_values, anchor_index)
            self._printDustConfiguration()

            if dry_run:
                self._print_dry_run(reader, labelled, time_values, anchor_index)
                return []

            for record in labelled:
                print(f"[macda2w] processing time_index={record.time_index}")
                # 每次处理一个时刻；start_index 是固定的日期标签和 xfcst 起算锚点。
                out_path = self._write_one_time(reader, record, time_values,
                                                anchor_index, lat, lon, sigma)
                written.append(out_path)
        # 离开 with 后输入文件关闭；written 仅包含已经写出并回读检查通过的路径。
        return written

    def _anchor_index(self, ntimes: int) -> int:
        """Return the stable label and XFCST anchor, which is ``start_index``.

        The anchor never depends on which time indices a single command selects,
        so a partial rerun writes exactly the filenames a full run writes.
        """

        anchor_index = self.cfg.start_index
        if anchor_index < 0 or anchor_index >= ntimes:
            raise ValueError(f"start_index {anchor_index} outside available time range "
                             f"0..{ntimes - 1}; it is the label and XFCST anchor")
        return anchor_index

    def _select_indices(self, ntimes: int, anchor_index: int,
                        time_indices: Iterable[int] | None,
                        max_times_override: int | None) -> list[int]:
        if time_indices is None:
            apply_config_limit = max_times_override is None
            indices = self._time_indices(ntimes, apply_config_limit=apply_config_limit)
        else:
            requested = [int(index) for index in time_indices]
            indices = sorted(dict.fromkeys(requested))
            if indices != requested:
                print(f"[macda2w] explicit time indices sorted and deduplicated: "
                      f"{indices}")
        if max_times_override is not None:
            if max_times_override < 1:
                raise ValueError("max_times must be at least 1")
            indices = indices[:max_times_override]
        self._validate_indices(indices, ntimes, anchor_index)
        return indices

    @staticmethod
    def _validate_indices(indices: list[int], ntimes: int, anchor_index: int) -> None:
        # 空列表、越界索引和早于锚点的索引均报错；排序去重由 _select_indices 负责。
        if not indices:
            raise ValueError("No MACDA time indices selected")
        invalid = [index for index in indices if index < 0 or index >= ntimes]
        if invalid:
            raise ValueError(f"time indices {invalid} outside available range "
                             f"0..{ntimes - 1}")
        early = [index for index in indices if index < anchor_index]
        if early:
            raise ValueError(f"time indices {early} come before the anchor record "
                             f"start_index={anchor_index}, so XFCST would be negative; "
                             f"lower start_index in the config to move the anchor")

    def _time_indices(self, ntimes: int, apply_config_limit: bool = True) -> list[int]:
        # end_index 为空时选到末帧；过大的终点截到可用末帧，起点不合法则直接拒绝。
        start = self.cfg.start_index
        end = self.cfg.end_index if self.cfg.end_index is not None else ntimes - 1
        end = min(end, ntimes - 1)
        if start < 0 or start >= ntimes:
            raise ValueError(f"start_index {start} outside available time range "
                             f"0..{ntimes - 1}")
        if end < start:
            raise ValueError(f"end_index {end} is before start_index {start}")
        # Python range 右端不含，因此加 1 以实现配置所声明的闭区间。
        indices = list(range(start, end + 1))
        if apply_config_limit and self.cfg.max_times is not None:
            indices = indices[: self.cfg.max_times]
        return indices

    def _build_labels(self, reader: MacdaReader, indices: list[int],
                      time_values: np.ndarray,
                      anchor_index: int) -> list[LabelledRecord]:
        """Label the whole selection against the anchor and validate the sequence."""

        anchor_date = reader.read_time_string(anchor_index)
        records = [MacdaRecord(time_index=index,
                               mars_date=reader.read_time_string(index),
                               time_sols=float(time_values[index]))
                   for index in indices]
        tolerance = self.cfg.time_tolerance_seconds
        labelled = build_label_sequence(records, strategy=self.cfg.hdate_strategy,
                                        anchor_date=anchor_date,
                                        tolerance_seconds=tolerance)
        short_years = short_mars_years_left(labelled)
        if short_years and self.cfg.hdate_strategy not in ANCHORED_HDATE_STRATEGIES:
            print(f"[macda2w] warning: this selection leaves MY {short_years}, which is "
                  f"shorter than the fixed {MARSWRF_FIXED_YEAR_SOLS}-sol MarsWRF year; "
                  f"metgrid will ask for a sol these labels do not contain")
        return labelled

    def _write_one_time(self, reader: MacdaReader, record: LabelledRecord,
                        time_values: np.ndarray, anchor_index: int, lat: np.ndarray,
                        lon: np.ndarray, sigma: np.ndarray) -> Path:
        time_index = record.time_index
        # 日期标签由完整序列统一生成；火星日期不能用地球 datetime 代替。
        hdate = record.hdate
        stamp = record.filename_stamp
        if record.relabelled:
            print(f"[macda2w] {self.cfg.hdate_strategy} relabel "
                  f"{record.mars_date.filename_stamp} -> {stamp}")
        if self.cfg.validate_time_alignment:
            # 可选核验字符串日期与数值 sol 时间是否一致，允许误差由配置的秒容差决定。
            validate_macda_time(record.mars_date, float(time_values[time_index]),
                                tolerance_seconds=self.cfg.time_tolerance_seconds)
        output_path = self.cfg.output_root / f"{self.cfg.output_prefix}:{stamp}"
        xfcst = self._xfcst_hours(time_values, time_index, anchor_index)

        # psfc 保留在排序后的源网格，用于 sigma 压力和柱光学厚度归一化，不能提前插到目标网格。
        psfc_raw = reader.read("psurf", time_index=time_index)
        psfc, src_lats, src_lons = orient_lat_lon(psfc_raw, lat, lon,
                                                  self.cfg.lon_convention)
        # 此处只需 psfc；每个物理字段还会取回其配套的已排序经纬度。
        del src_lats, src_lons
        self._printDustOpacityDiagnostic(reader, time_index, sigma, psfc_raw)
        expectedDustRecords = self._expectedDustRecords(reader, time_index, lat,
                                                         lon, sigma, psfc)

        # 一个时刻对应一个中间文件；writer 在文件中顺序写入各字段、各压力层记录。
        with WrfIntermediateWriter(output_path, self.grid, hdate=hdate,
                                   map_source=self.cfg.map_source,
                                   xfcst=xfcst) as writer:
            for field in self.fields:
                # 禁用行始终跳过；可选行还要由全局 emit_optional_fields 允许才输出。
                if not field.enabled:
                    continue
                if not field.required and not self.cfg.emit_optional_fields:
                    continue
                self._write_field(reader, writer, field, time_index,
                    lat, lon, sigma, psfc)
        # 写入器关闭后回读，核对日期、网格尺寸、预报时效以及全部启用的必需字段。
        required_fields = {field.aim_v for field in self.fields
                           if field.enabled and field.required}
        records = validate_intermediate_file(output_path, expected_hdate=hdate,
            expected_shape=(self.grid.nlat, self.grid.nlon), expected_xfcst=xfcst,
            required_fields=required_fields)
        self._validateDustRecords(records, expectedDustRecords)
        return output_path

    def _write_field(self, reader: MacdaReader, writer: WrfIntermediateWriter,
                     field: FieldSpec, time_index: int, lat: np.ndarray,
                     lon: np.ndarray, sigma: np.ndarray, psfc: np.ndarray) -> None:
        # 人工压力场不从源变量读数：每个目标层的二维 slab 都填入该层压力，单位 Pa。
        if field.kind == "3d_pressure":
            for plev in self.cfg.plev_pa:
                slab = np.full((self.grid.nlat, self.grid.nlon), plev,
                               dtype=np.float32)
                writer.write_field(field.aim_v, slab, field.units, field.desc,
                                   xlvl=plev)
            return

        if field.kind == "3d_constant":
            # 当前常量分支固定输出零，不从变量表读取常数值，也不能解释为反演出的尘埃质量。
            for plev in self.cfg.plev_pa:
                slab = np.zeros((self.grid.nlat, self.grid.nlon), dtype=np.float32)
                writer.write_field(field.aim_v, slab, field.units, field.desc,
                                   xlvl=plev)
            return

        if not reader.has_var(field.src_v):
            # 必需源变量缺失就停止；可选变量缺失只打印跳过提示，不伪造填充值。
            if field.required:
                raise KeyError(f"Required MACDA variable missing: {field.src_v}")
            print(f"[macda2w] skip optional missing variable {field.src_v}")
            return

        if field.kind == "2d":
            # 二维场先取单时刻，再统一方向，使其与上面的源网格 psfc 逐格对应。
            raw = reader.read(field.src_v, time_index=time_index)
            slab, src_lats, src_lons = orient_lat_lon(raw,
                lat, lon, self.cfg.lon_convention)
            # 在源网格完成物理变换后才水平插值；否则 psfc 与光学厚度可能不在同一网格。
            slab = self._transform_2d(field, slab, psfc)
            out = horizontal_interp_2d(slab, src_lats, src_lons, self.grid)
            self._ensure_finite(field, out)
            # 地表字段使用特殊层码 XLVL_SURFACE；UU/VV 声明为地理东/北方向的风。
            writer.write_field(field.aim_v, out, field.units, field.desc,
                xlvl=XLVL_SURFACE,
                is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}),)
            return

        if field.kind == "3d_sigma":
            # 三维源场为 (sigma, lat, lon)，与 psfc 使用同样的经纬度排序。
            raw = reader.read(field.src_v, time_index=time_index)
            volume, src_lats, src_lons = orient_lat_lon(raw, lat, lon,
                                                        self.cfg.lon_convention)
            volume = self._transform_3d(field, volume)
            # MACDA 的 p=sigma*ps 由垂直模块处理；此处不能直接代入 MCD 的 A+B*ps 混合层。
            pvol = interp_sigma_to_pressure(volume, psfc, sigma, self.cfg.plev_pa)
            for idx, plev in enumerate(self.cfg.plev_pa):
                # 每个压力层独立水平插值，再检查有限性，最后写出一条二维 WPS 记录。
                out = horizontal_interp_2d(pvol[idx], src_lats, src_lons, self.grid)
                self._ensure_finite(field, out)
                writer.write_field(field.aim_v, out, field.units, field.desc,
                                   xlvl=plev,
                                   is_wind_earth_rel=int(field.aim_v in {"UU", "VV"}))
            return

        # CSV 解析器允许 3d_hybrid 供 MCD 复用，但这个 MACDA 写场方法未实现其计算分支。
        # 不认识的 kind 明确报错，避免默默输出未经正确垂直转换的数据。
        raise ValueError(f"Unsupported field kind {field.kind} for {field.aim_v}")

    @staticmethod
    def _ensure_finite(field: FieldSpec, data: np.ndarray) -> None:
        # 统计 NaN 与正负无穷的总数；拒绝带无效值的场进入中间文件，不用零或近邻遮盖。
        bad = int(np.size(data) - np.isfinite(data).sum())
        if bad:
            raise ValueError(f"Field {field.aim_v} contains {bad} non-finite values "
                             f"after conversion")
        if field.transform in DUST_TRANSFORMS and np.any(data < 0.0):
            minimum = float(np.min(data))
            raise ValueError(f"Field {field.aim_v} contains negative dust mixing "
                             f"ratios; minimum={minimum:.8g}")

    def _transform_2d(self, field: FieldSpec, slab: np.ndarray,
                      psfc: np.ndarray) -> np.ndarray:
        # identity 表示源数值已符合目标定义；变量表本身不会根据 units 自动做量纲转换。
        if field.transform == "identity":
            return slab
        if field.transform == "dust_to_tau_7mb":
            # 名称保留旧接口的 7mb，但实际参考压力来自 cfg.tau_reference_pressure_pa。
            # 数学关系为 tau_ref=tau_column*p_ref/ps；它不改变热力基态 P0=610 Pa。
            # ps<=0 输出 NaN，再由有限性检查拒绝；errstate 只关闭警告，不补造有效值。
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(psfc > 0.0,
                                slab / psfc * self.cfg.tau_reference_pressure_pa,
                                np.nan)
        raise ValueError(f"Unsupported 2-D transform {field.transform}")

    def _transform_3d(self, field: FieldSpec, volume: np.ndarray) -> np.ndarray:
        if field.transform == "identity":
            return volume
        if field.transform == "geop_to_height":
            # 位势 Phi 的单位 m²/s²，除以火星重力得到位势高度 GHT，单位 m。
            return volume / self.cfg.mars_gravity
        if field.transform == DUST_BIN1_TRANSFORM:
            return tracerFromDustMass(volume, self.cfg.dustMassFraction1,
                                      self.cfg.dustWeight1)
        if field.transform == DUST_BIN2_TRANSFORM:
            return tracerFromDustMass(volume, self.cfg.dustMassFraction2,
                                      self.cfg.dustWeight2)
        raise ValueError(f"Unsupported 3-D transform {field.transform}")

    def _expectedDustRecords(self, reader: MacdaReader, timeIndex: int,
                             lat: np.ndarray, lon: np.ndarray, sigma: np.ndarray,
                             surfacePressure: np.ndarray) -> list[tuple[float, np.ndarray]]:
        dustFields = [field for field in self.fields
                      if field.enabled and field.transform in DUST_TRANSFORMS]
        if not dustFields:
            return []

        sourceVariable = dustFields[0].src_v
        raw = reader.read(sourceVariable, time_index=timeIndex)
        volume, srcLats, srcLons = orient_lat_lon(raw, lat, lon, self.cfg.lon_convention)
        pressureVolume = interp_sigma_to_pressure(volume, surfacePressure, sigma, self.cfg.plev_pa)
        expectedRecords: list[tuple[float, np.ndarray]] = []
        for index, pressureLevel in enumerate(self.cfg.plev_pa):
            slab = horizontal_interp_2d(pressureVolume[index], srcLats, srcLons, self.grid)
            self._ensure_finite(dustFields[0], slab)
            expectedRecords.append((pressureLevel, slab))
        return expectedRecords

    def _validateDustRecords(self, records: list[IntermediateRecord],
                              expectedDustRecords: list[tuple[float, np.ndarray]] | None = None) -> None:
        bin1List = [record for record in records if record.field == "TRC01"]
        bin2List = [record for record in records if record.field == "TRC02"]
        if not bin1List and not bin2List:
            return
        expectedCount = len(self.cfg.plev_pa)
        if len(bin1List) != expectedCount or len(bin2List) != expectedCount:
            raise ValueError(f"TRC01/TRC02 must each contain {expectedCount} pressure levels; "
                             f"found {len(bin1List)}/{len(bin2List)}")

        bin1Records = self._sortDustRecords(bin1List, "TRC01")
        bin2Records = self._sortDustRecords(bin2List, "TRC02")
        bin1Levels = np.asarray([record.xlvl for record in bin1Records], dtype=np.float64)
        bin2Levels = np.asarray([record.xlvl for record in bin2Records], dtype=np.float64)
        if not np.allclose(bin1Levels, bin2Levels, rtol=1.0e-6, atol=1.0e-6):
            raise ValueError("TRC01 and TRC02 do not contain the same pressure levels")
        expectedLevels = np.asarray(self.cfg.plev_pa, dtype=np.float64)
        if not np.allclose(np.sort(expectedLevels), bin1Levels, rtol=1.0e-6, atol=1.0e-6):
            raise ValueError("TRC01/TRC02 pressure levels do not match plev_pa")
        if expectedDustRecords is not None:
            if len(expectedDustRecords) != expectedCount:
                raise ValueError("Expected dust records do not match the configured pressure levels")
            expectedDustRecords = sorted(expectedDustRecords, key=lambda record: record[0])
            expectedLevels = np.asarray([record[0] for record in expectedDustRecords],
                                        dtype=np.float64)
            if not np.allclose(expectedLevels, bin1Levels, rtol=1.0e-6, atol=1.0e-6):
                raise ValueError("Expected dust records do not match TRC01/TRC02 pressure levels")

        combinedSlabs: list[np.ndarray] = []
        for index, (bin1Record, bin2Record) in enumerate(zip(bin1Records, bin2Records)):
            level = bin1Record.xlvl
            tracer1 = bin1Record.slab
            tracer2 = bin2Record.slab
            if np.any(tracer1 < 0.0) or np.any(tracer2 < 0.0):
                raise ValueError(f"TRC01/TRC02 contain negative values at {level:.8g} Pa")
            combined = combineDustTracers(tracer1, tracer2, self.cfg.dustWeight1,
                self.cfg.dustWeight2)
            combinedSlabs.append(combined)
            sourceEstimates: list[np.ndarray] = []
            if self.cfg.dustMassFraction1 > 0.0:
                sourceEstimates.append(tracer1 * self.cfg.dustWeight1 / self.cfg.dustMassFraction1)
            if self.cfg.dustMassFraction2 > 0.0:
                sourceEstimates.append(tracer2 * self.cfg.dustWeight2 / self.cfg.dustMassFraction2)
            for sourceEstimate in sourceEstimates:
                if not np.allclose(combined, sourceEstimate, rtol=1.0e-6, atol=1.0e-15):
                    raise ValueError(f"TRC01/TRC02 do not conserve dust mass at {level:.8g} Pa")
            if expectedDustRecords is not None:
                sourceDust = expectedDustRecords[index][1]
                try:
                    validateDustPartition(sourceDust, tracer1, tracer2,
                                          self.cfg.dustWeight1, self.cfg.dustWeight2)
                except ValueError as error:
                    raise ValueError(f"TRC01/TRC02 do not reconstruct source dust at {level:.8g} Pa") from error

        combinedVolume = np.stack(combinedSlabs)
        if not np.any(combinedVolume > 0.0):
            raise ValueError("TRC01/TRC02 are zero at every pressure level")
        print(f"[macda2w] dust mass verified: levels={len(combinedSlabs)} "
              f"min={float(np.min(combinedVolume)):.8g} "
              f"max={float(np.max(combinedVolume)):.8g}")

    @staticmethod
    def _sortDustRecords(records: list[IntermediateRecord], fieldName: str) -> list[IntermediateRecord]:
        sortedRecords = sorted(records, key=lambda record: record.xlvl)
        levels = np.asarray([record.xlvl for record in sortedRecords], dtype=np.float64)
        if np.any(np.isclose(np.diff(levels), 0.0, rtol=1.0e-6, atol=1.0e-6)):
            raise ValueError(f"{fieldName} contains duplicate pressure levels")
        return sortedRecords

    def _printDustConfiguration(self) -> None:
        if not any(field.transform in DUST_TRANSFORMS for field in self.fields):
            return
        scale = opticalDepthScaling(self.cfg.macdaDustDensity, self.cfg.macdaDustEffectiveRadius,
            self.cfg.macdaDustExtinctionEfficiency, self.cfg.mars_gravity,
            self.cfg.marswrfRadiationDustDensity, self.cfg.marswrfRadiationGravity,
            self.cfg.marswrfDustDiameter1, self.cfg.marswrfDustDiameter2,
            self.cfg.marswrfDustExtinctionEfficiency1, self.cfg.marswrfDustExtinctionEfficiency2,
            self.cfg.dustMassFraction1, self.cfg.dustMassFraction2)
        print(f"[macda2w] dust partition: mass fractions="
              f"{self.cfg.dustMassFraction1:.6g}/{self.cfg.dustMassFraction2:.6g}, "
              f"MarsWRF wt={self.cfg.dustWeight1:.6g}/{self.cfg.dustWeight2:.6g}")
        print(f"[macda2w] estimated MarsWRF/MACDA dust-opacity scaling={scale:.6g}")
        if not np.isclose(scale, 1.0, rtol=0.05, atol=0.0):
            print("[macda2w] warning: mass is conserved, but the configured MarsWRF "
                  "particle optics cannot reproduce MACDA opacity within 5%")

    def _printDustOpacityDiagnostic(self, reader: MacdaReader, timeIndex: int, sigma: np.ndarray,
        surfacePressure: np.ndarray) -> None:
        if not any(field.transform in DUST_TRANSFORMS for field in self.fields):
            return
        if not reader.has_var("dustmmr") or not reader.has_var("coldust"):
            return

        dust = reader.read("dustmmr", time_index=timeIndex)
        observed = reader.read("coldust", time_index=timeIndex)
        macdaMassExtinction = massExtinctionCoefficient(self.cfg.macdaDustExtinctionEfficiency,
            self.cfg.macdaDustDensity, self.cfg.macdaDustEffectiveRadius)
        macdaOpticalDepth = columnOpticalDepth(dust, surfacePressure, sigma, macdaMassExtinction,
            self.cfg.mars_gravity)
        scale = opticalDepthScaling(self.cfg.macdaDustDensity, self.cfg.macdaDustEffectiveRadius,
            self.cfg.macdaDustExtinctionEfficiency, self.cfg.mars_gravity,
            self.cfg.marswrfRadiationDustDensity, self.cfg.marswrfRadiationGravity,
            self.cfg.marswrfDustDiameter1, self.cfg.marswrfDustDiameter2,
            self.cfg.marswrfDustExtinctionEfficiency1, self.cfg.marswrfDustExtinctionEfficiency2,
            self.cfg.dustMassFraction1, self.cfg.dustMassFraction2)
        valid = np.isfinite(observed) & np.isfinite(macdaOpticalDepth) & (observed > 0.0)
        if not np.any(valid):
            print("[macda2w] warning: coldust has no positive finite values for closure")
            return
        macdaRatio = float(np.median(macdaOpticalDepth[valid] / observed[valid]))
        marswrfRatio = macdaRatio * scale
        print(f"[macda2w] opacity closure: MACDA-derived/coldust median={macdaRatio:.6g}, "
              f"MarsWRF-configured/coldust median={marswrfRatio:.6g}")

    def _xfcst_hours(self, time_values: np.ndarray, time_index: int,
                     anchor_index: int) -> float:
        """Forecast hours since the anchor record, in Martian hours."""

        # xfcst 仅写入文件头；关掉起算时效时所有记录都标 0。
        if not self.cfg.xfcst_from_start:
            return 0.0
        # MACDA 数值 time 以 sol 计，乘 24 得火星模型标签小时，不是地球小时。
        hours = float((time_values[time_index] - time_values[anchor_index]) * 24.0)
        if hours < 0.0:
            raise ValueError(f"XFCST for time index {time_index} is {hours:.6g} h, "
                             f"before the anchor record start_index={anchor_index}")
        return hours

    def _print_dry_run(self, reader: MacdaReader, labelled: list[LabelledRecord],
                       time_values: np.ndarray, anchor_index: int) -> None:
        # 打印输入后端、网格、压力层，以及所选每个时刻的日期标签和预报时效。
        anchor_date = reader.read_time_string(anchor_index)
        relabelled = [record for record in labelled if record.relabelled]
        print(f"input_file={self.cfg.input_file}")
        print(f"reader_backend={reader.backend}")
        print(f"output_root={self.cfg.output_root}")
        print(f"target_grid={self.grid.nlat}x{self.grid.nlon}")
        print(f"plev_pa={self.cfg.plev_pa}")
        print(f"hdate_strategy={self.cfg.hdate_strategy}")
        print(f"label_anchor=time[{anchor_index}] {anchor_date}")
        print(f"records={len(labelled)} relabelled={len(relabelled)} "
              f"first={labelled[0].hdate[:19]} last={labelled[-1].hdate[:19]}")
        print("times:")
        for record in labelled:
            if self.cfg.validate_time_alignment:
                validate_macda_time(record.mars_date, record.time_sols,
                                    tolerance_seconds=self.cfg.time_tolerance_seconds)
            xfcst = self._xfcst_hours(time_values, record.time_index, anchor_index)
            mark = " relabelled" if record.relabelled else ""
            print(f"  time[{record.time_index}]={record.time_sols:.12g} sol "
                  f"{record.mars_date.macda_string} -> {record.mars_date.wrf_date} "
                  f"-> label {record.hdate[:19]} xfcst={xfcst:.6g} h{mark}")
        print("fields:")
        for field in self.fields:
            # 计划显示的字段使用与实际写出完全相同的启用/必需/可选筛选条件。
            if field.enabled and (field.required or self.cfg.emit_optional_fields):
                status = "required" if field.required else "optional"
                print(f"  {field.aim_v:10s} <- {field.src_v:12s} "
                      f"{field.kind:13s} {status}")


def read_field_table(path: str | Path) -> list[FieldSpec]:
    # 共用 CSV 读取入口：第一行为列名，每一数据行建立不可变 FieldSpec，按原行顺序保留。
    fields: list[FieldSpec] = []
    with Path(path).open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            # 清除字符串两端空格；布尔列接受的写法由 _csv_bool 统一，避免各入口自行解析。
            field = FieldSpec(src_v=row["src_v"].strip(), aim_v=row["aim_v"].strip(),
                              units=row["units"].strip(), kind=row["kind"].strip(),
                              transform=row["transform"].strip(),
                              required=_csv_bool(row["required"]),
                              enabled=_csv_bool(row["enabled"]),
                              desc=row["desc"].strip())
            # 源/目标名不可为空；WPS 固定长度字段名最多 9 个字符。
            if not field.src_v or not field.aim_v:
                raise ValueError(f"Empty src_v/aim_v in variable table {path}")
            if len(field.aim_v) > 9:
                raise ValueError(f"WRF intermediate field name exceeds 9 characters: "
                                 f"{field.aim_v}")
            # 允许 MCD 的 3d_hybrid 被共用读取器接受，但执行它的是 MCD pipeline，不是本类。
            if field.kind not in {"2d", "3d_sigma", "3d_pressure", "3d_constant", "3d_hybrid"}:
                raise ValueError(f"Unsupported kind {field.kind!r} for {field.aim_v}")
            fields.append(field)
    # 空表没有任何可转换字段，直接拒绝，不生成内容为空但看似成功的中间文件。
    if not fields:
        raise ValueError(f"Variable table contains no fields: {path}")
    return fields


def validateDustFieldTable(fields: list[FieldSpec]) -> None:
    """Require a complete, consistently configured two-bin MACDA dust mapping."""

    dustFields = [field for field in fields
                  if field.enabled and field.transform in DUST_TRANSFORMS]
    if not dustFields:
        return
    expectedTransforms = {DUST_BIN1_TRANSFORM, DUST_BIN2_TRANSFORM}
    transforms = {field.transform for field in dustFields}
    if transforms != expectedTransforms or len(dustFields) != 2:
        raise ValueError("Dust field table must contain exactly one mapping for each dust bin")
    if {field.aim_v for field in dustFields} != {"TRC01", "TRC02"}:
        raise ValueError("Dust field table must map the two dust bins to TRC01 and TRC02")
    if any(field.kind != "3d_sigma" for field in dustFields):
        raise ValueError("Dust field table must convert dust bins from sigma levels")
    if any(not field.required for field in dustFields):
        raise ValueError("Dust field table must require both dust bins")
    if len({field.src_v for field in dustFields}) != 1:
        raise ValueError("Dust field table must use one source variable for both dust bins")


def _csv_bool(value: str) -> bool:
    # 去空格、转小写后匹配肯定集合；不在集合中的值都为 False，原有解析行为保持不变。
    return value.strip().lower() in {"1", "true", "yes", "on"}
