"""Continuous MCD background preparation; not an observation assimilation system.

The MACDA grid/orientation, field table, calendar and binary WPS writer are
reused. MCD-specific season/solar sampling and hybrid/dust physics are separate.
No input field here proves delivery into real.exe or a WRF radiation scheme.
"""

# 本模块把 MCD 读取器、时间坐标、垂直坐标和沙尘诊断组合成可检查的背景产品。
# build_state 构造一个时刻的场；write_states 将多个时刻逐个写入 NetCDF。
# export_intermediate 是额外的 WPS 输出入口，不等同于已经生成或验证 WRF 初始场。
# 月平均插值、可选合成 EOF 扰动均不构成观测同化，因此输出明确标记为非再分析。

# asdict/JSON 用于把物理参数完整记入文件；contextmanager 与 UUID 用于安全发布输出。
from dataclasses import asdict
from contextlib import contextmanager
import json
import os
from pathlib import Path
from uuid import uuid4

# Dataset 管理 NetCDF 文件和属性；NumPy 处理各层、各格点的数组运算。
from netCDF4 import Dataset
import numpy as np

# 复用 MACDA 的变量表、水平坐标排序和 WPS 二进制写出功能。
from .converter import read_field_table
from .grid import RegularLatLonGrid, orient_lat_lon
# MCD 专用模块负责混合压力、静力高度、尘埃分 bin、火星时序和显式扰动。
from .mcd_coordinate import hybrid_pressure, pressure_interfaces, hydrostatic_height
from .mcd_coordinate import interpolate_to_pressure
from .mcd_dust import DustBinParameters, column_dust_diagnostics
from .mcd_time import iter_mcd_times, SI_SECONDS_PER_SOL, MCD_ORBIT_SOLS
from .mcd_reader import McdReader
from .mcd_value import apply_anomaly
from .wrf_intermediate import WrfIntermediateWriter, validate_intermediate_file


# 只改输出字段名称，不在这里改数值；源 tau_pref_gcm 明确标为可见光、610 Pa 光学厚度。
SOURCE_NAMES = {"temp": "temperature", "tau_pref_gcm": "tau_visible610_mcd"}
# NetCDF 输出单位表。单位字符串统一后，后续独立检查不必再猜测变量的物理定义。
UNITS = {"temperature": "K", "u": "m s-1", "v": "m s-1", "rho": "kg m-3",
    "dustq": "kg kg-1", "reffdust": "m", "vmr_h2ovapor": "mol mol-1", "ps": "Pa",
    "tsurf": "K", "co2ice": "kg m-2", "h2oice": "kg m-2", "pressure": "Pa",
    "pressure_interface": "Pa", "ght": "m", "height_interface": "m", "eta_diagnostic": "1",
    "q_bin1": "kg kg-1", "q_bin2": "kg kg-1", "reff_matched": "m", "reff_residual": "m",
    "radius_clipped": "1", "radiation_mass_residual": "kg kg-1", "specific_humidity": "kg kg-1",
    "effective_gas_constant": "J kg-1 K-1", "rms_temperature": "K", "rms_u": "m s-1",
    "rms_v": "m s-1", "rms_ps": "Pa", "rms_tau_visible610_mcd": "1"}
# 只从底层沙尘诊断中挑选这些字段输出：物理质量、粒径残差、柱质量和不同光学假设。
# q_bin1/q_bin2 是物理质量分数；辐射权重或下限产生的变化另存 residual，不反改物理质量。
DUST_FIELDS = ("q_bin1", "q_bin2", "reff_matched", "reff_residual", "radius_clipped",
    "radiation_mass_residual", "column_dust_mass", "column_dust_mass_bin1", "column_dust_mass_bin2",
    "tau_visible_spherical", "tau_visible610_spherical", "tau_ir_abs610_assumed",
    "tau_bins_visible610_unweighted", "tau_wrf_native", "tau_wrf_native610", "tau_wrf_native700",
    "tau_wrf_weighting_residual610")
# 仅允许这个量级的地表冰负舍入尾数归零；真实负冰量必须拒绝，不能整体截断掩盖。
ICE_ROUNDOFF_TOLERANCE_KG_M2 = 1.0e-18


@contextmanager
def new_dataset(path):
    """Publish only closed complete products; never overwrite an existing file."""
    # 统一路径并创建父目录；已有正式文件一律拒绝覆盖。
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    # 临时文件与目标位于同一目录，UUID 避免多个任务使用相同临时名称。
    temporary = path.with_name(f"{path.name}.{uuid4().hex}.partial")
    try:
        with Dataset(temporary, "w", format="NETCDF4") as output:
            # yield 把可写对象交给 with 调用方；调用方完成后先关闭文件，再发布正式名称。
            yield output
        # 硬链接让完整内容以正式名称出现；即使竞争者刚创建目标，此操作也不会覆盖它。
        os.link(temporary, path)
    finally:
        # 成功时移除多余临时名称，异常时删除不完整临时文件；正式目标不在清理范围内。
        temporary.unlink(missing_ok=True)


def build_state(reader, point, dust_parameters=None, variability=None):
    """Return oriented native-grid state and transparent physical diagnostics."""
    # source 是指定时刻的插值场字典；大气三维量为 (level, lat, lon)，地表量为 (lat, lon)。
    source = reader.sample(point)
    ice_corrections = {}
    for name in ("co2ice", "h2oice"):
        values = source[name]
        # 超过明确舍入容差的负值不允许进入背景场，也不把它当成可接受的冰量。
        if np.any(values < -ICE_ROUNDOFF_TOLERANCE_KG_M2):
            raise ValueError(f"{name}: negative ice exceeds explicit roundoff tolerance")
        # correction 保存所加的正修正量，可由 original=输出冰量-correction 恢复原数组。
        ice_corrections[f"{name}_roundoff_correction"] = np.maximum(-values, 0.0)
        source[name] = np.maximum(values, 0.0)
    # RMS、冰量和光学厚度须有限且非负；RMS 只表示变化尺度，不在此充当距平扣除。
    for name, array in source.items():
        if name.startswith("rms_") or name in ("co2ice", "h2oice", "tau_pref_gcm"):
            if np.any(array < 0) or not np.all(np.isfinite(array)):
                raise ValueError(f"{name}: expected finite nonnegative source values")
    if np.any(source["vmr_h2ovapor"] < 0) or np.any(source["vmr_h2ovapor"] > 1):
        raise ValueError("water vapour mole fraction must lie within [0, 1]")
    # p=A+B*ps 还原各柱层中心的实际压力；用同时刻 p、rho、T 推算源背景的有效气体常数。
    # 这个 R 保留插值背景的状态方程关系，不能把地球空气的气体常数直接代入火星。
    background_pressure = hybrid_pressure(source["ps"], reader.aps, reader.bps)
    background_r = background_pressure / (source["rho"] * source["temp"])
    increments = {}
    if variability is not None:
        # 可选 EOF 提供有符号且有空间结构的合成扰动；它不是由 RMS 相减得到的观测真值。
        perturbation = variability.sample(point.elapsed_sols, reader.latitude,
            reader.longitude, reader.aps, reader.bps, source_ps=source["ps"])
        for name in ("temp", "u", "v", "ps"):
            # 单独保留增量，以便区分原背景与候选场；apply_anomaly 按声明的 signed_anomaly 相加。
            increments[f"increment_{SOURCE_NAMES.get(name, name)}"] = perturbation[name]
            source[name] = apply_anomaly(source[name], perturbation[name], "signed_anomaly")
    # 地表压力可能已被扰动，因此必须重新计算整根混合压力柱。
    pressure = hybrid_pressure(source["ps"], reader.aps, reader.bps)
    if np.any(source["temp"] <= 0) or np.any(source["rho"] <= 0):
        raise ValueError("nonphysical temperature or density; candidate rejected, not clipped")
    if variability is not None:
        # 扰动后用原有效 R 重算密度，保持 p=rho*R*T；不独立给密度添加另一个不一致扰动。
        source["rho"] = pressure / (background_r * source["temp"])
    # 层界面缺少原生边界时由中心压力推断：底界为 ps，内部为对数压力中点，顶界为正值外推。
    # 这些界面给出选定离散柱的 dp；不能宣称它们就是 PCM 原生有限体积层的真实边界。
    interface = pressure_interfaces(pressure, source["ps"])
    # 从源地形向上按分层等温静力关系积分高度。areoid_radius 是半径，orography 才是地形高。
    # density 参数使高度计算采用 p/(rho*T) 的有效 R，输出 GHT 以当地 areoid 为高度基准。
    height = hydrostatic_height(pressure, interface.pressure_pa, source["temp"],
        reader.orography, reader.areoid_radius, density_kg_m3=source["rho"])
    # 汇总原变量、微小冰量修正和合成增量；只在这里统一 temp 等输出名称。
    fields = {SOURCE_NAMES.get(name, name): value for name, value in source.items()}
    fields.update(ice_corrections)
    fields.update(increments)
    # eta_diagnostic=(p-ptop)/(ps-ptop) 仅是总压诊断坐标。
    # MarsWRF 真正的干静力 eta、目标地形和模式层仍由接收链及 real.exe 构建。
    fields.update({"pressure": pressure, "pressure_interface": interface.pressure_pa,
        "ght": height.geopotential_height_areoid_m,
        "height_interface": height.interface_geopotential_height_areoid_m,
        "effective_gas_constant": height.effective_gas_constant,
        "eta_diagnostic": (pressure - interface.pressure_pa[-1]) / (source["ps"] - interface.pressure_pa[-1])})
    # 逐层质量分数乘 dp/g 得柱质量；用源有效半径分配两个物理 bin，并计算不同光学诊断。
    # 输入 dustq 已是 kg/kg，粒径只影响分 bin 与光学计算，不凭粒径重新“创造”质量。
    dust = column_dust_diagnostics(source["dustq"], source["reffdust"], interface.pressure_pa, dust_parameters)
    fields.update({name: dust[name] for name in DUST_FIELDS})
    # With MCD effective R, mean molecular mass is R_universal/R_effective.
    # VMR * Mwater/Mair is mass fraction; unlike dry mixing ratio this is q/(1+q).
    # 水汽摩尔分数乘 Mwater/Mair 转为水质量/总空气质量；Mair=Ru/R，有效 R 保留火星气体组成。
    # 0.01801528 为水的 kg/mol 摩尔质量，8.314462618 为通用气体常数；结果是比湿而非干混合比。
    fields["specific_humidity"] = source["vmr_h2ovapor"] * 0.01801528 * height.effective_gas_constant / 8.314462618
    if np.any(fields["specific_humidity"] < 0) or np.any(fields["specific_humidity"] >= 1):
        raise ValueError("invalid water mass fraction")
    for name, array in fields.items():
        valid = np.isfinite(array)
        if name in ("reff_matched", "reff_residual"):
            # 无尘处的有效半径没有物理定义，只允许这些半径诊断在 dustq=0 处保留 NaN。
            valid = valid | ((source["dustq"] == 0) & np.isnan(array))
        if not np.all(valid):
            raise ValueError(f"{name}: missing/nonfinite physical background; sample rejected")
        # 对每个字段执行相同的纬度升序、经度排序，保证层中心、界面与地表的水平索引一致。
        fields[name] = orient_lat_lon(array, reader.latitude, reader.longitude, "native")[0]
    # 保存重建方法以区分源量和推断量；EOF 设置及相位也作为可追溯信息返回。
    metadata = {"interface_method": interface.method, "top_boundary": interface.top_boundary,
        "height_method": height.method}
    if variability is not None:
        metadata["variability_metadata"] = json.dumps(perturbation["metadata"], sort_keys=True)
    return fields, metadata


def add_time_variables(dataset, points):
    # 一个 time 维承载多套互相可核对的时间表示，不能把火星标签秒误读为地球 UTC 秒。
    dataset.createDimension("time", len(points))
    # time 是 SI 经过秒；calendar_seconds 每 sol 采用 86400 个标签秒，仅服务 MarsWRF 日期。
    # model_hour 是本初子午线平太阳时，ltst0 是对应真太阳时，ls 表示季节而非均匀计时。
    # orbit_wrapped 标记超过 668.6-sol 理想轨道的末端记录，另一个标志标记月中心插值跨年。
    items = {"time": ("f8", [point.elapsed_seconds for point in points], "s"),
        "calendar_seconds": ("i8", [point.calendar_seconds for point in points], "model_label_seconds"),
        "sol_of_year": ("i4", [point.sol_of_year for point in points], "1"),
        "model_hour": ("f8", [point.mean_solar_time for point in points], "Mars_hour"),
        "ltst0": ("f8", [point.true_solar_time for point in points], "Mars_hour"),
        "ls": ("f8", [point.solar_longitude for point in points], "degree"),
        "orbit_wrapped": ("i1", [point.orbit_wrapped for point in points], "1"),
        "seasonal_boundary_cycled": ("i1", [point.seasonal_weights.crosses_year_boundary for point in points], "1")}
    for name, (dtype, data, units) in items.items():
        # 每个时间量独立写入类型、数组和单位，避免只靠日期字符串解释数值时间。
        variable = dataset.createVariable(name, dtype, ("time",))
        variable[:] = data
        variable.units = units
    dataset["time"].long_name = "SI elapsed seconds since MY start; not UTC or Earth CF calendar"
    dataset["calendar_seconds"].long_name = "86400 labelled seconds per sol for MarsWRF date strings"
    # wrf_date 保留 YYYY-DDDDD_HH:MM:SS 的火星年/年内 sol 标签，不套用地球日历。
    variable = dataset.createVariable("wrf_date", str, ("time",))
    variable[:] = np.asarray([point.wrf_date for point in points], dtype=object)


def global_metadata(dataset, reader, parameters, candidate=False):
    # 明确区分纯插值背景与合成 EOF 候选；两者都没有观测同化，也没有 WRF 接收证明。
    dataset.product_type = "mcd_synthetic_eof_candidate" if candidate else "mcd_interpolated_background"
    dataset.is_reanalysis = "false: no observation assimilation has been performed"
    dataset.mars_year = reader.year
    # 统一物理参考压力为 610 Pa，即 6.1 mb；独立光学诊断的历史 700 Pa 定义不改变 P0。
    dataset.reference_pressure_pa = 610.0
    dataset.gravity_m_s2 = parameters.gravity_m_s2
    dataset.time_standard = "prime-meridian LMST/MTC model calendar; zero epoch is assumed spring equinox"
    dataset.si_seconds_per_sol = SI_SECONDS_PER_SOL
    dataset.mcd_idealized_orbit_sols = MCD_ORBIT_SOLS
    # 跨年循环、总压 eta 和原生网格政策均写入属性，让后续使用者知道资料的重建边界。
    dataset.year_boundary_policy = "same-MY climatological cycling; not adjacent-year observations"
    dataset.eta_definition = "diagnostic total-pressure (p-ptop)/(ps-ptop); actual dry eta belongs to real.exe"
    dataset.source_directory = str(reader.source_dir)
    dataset.static_source = str(reader.mountain_file)
    dataset.horizontal_grid_policy = "native 64x49 for this archive; no resolution enhancement"
    # 参数整体序列化，包含粒径、密度、效率、权重、重力和下限，供续算及诊断核验。
    dataset.dust_parameters = json.dumps(asdict(parameters), sort_keys=True)
    dataset.rms_definition = "interpolated unsigned source variability scale; not signed anomaly or full B"
    dataset.optical_definition = "MCD visible 610-Pa tau; spherical and 9.3um assumptions labelled separately"
    dataset.wrf_delivery = "not verified: WPS/real/physics receipt is outside this preparation product"
    dataset.ice_roundoff_tolerance_kg_m2 = ICE_ROUNDOFF_TOLERANCE_KG_M2
    dataset.ice_roundoff_policy = "negative values within tolerance set to zero; raw=ice-roundoff_correction"


def write_calendar(path, reader, points, parameters):
    # 独立日历文件只写全局属性和时间坐标；空时间列表或已有目标文件均拒绝。
    if not points:
        raise ValueError("calendar requires at least one time point")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with new_dataset(path) as output:
        global_metadata(output, reader, parameters)
        add_time_variables(output, points)


def write_states(path, reader, points, parameters=None, variability=None):
    """Stream a block one time at a time; no full-year expanded RAM array."""
    # 调用者通常传入一个 sol 的时刻列表；这里只要求非空，不把全年所有三维场展开到内存。
    if not points:
        raise ValueError("state product requires at least one time point")
    parameters = parameters or DustBinParameters()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with new_dataset(path) as output:
        # 先建立产品级属性和时间坐标，再建立水平/垂直坐标，最后逐时写入物理场。
        global_metadata(output, reader, parameters, variability is not None)
        add_time_variables(output, points)
        for name, values in (("lat", np.sort(reader.latitude)), ("lon", np.sort(reader.longitude)),
            ("level", np.arange(len(reader.aps))), ("interface", np.arange(len(reader.aps) + 1))):
            # level 与 interface 是索引；真正压力随 ps 变化，另写三维 pressure 字段。
            output.createDimension(name, len(values))
            variable = output.createVariable(name, "f8", (name,))
            variable[:] = values
            variable.units = {"lat": "degrees_north", "lon": "degrees_east"}.get(name, "1")
        for name, values, units in (("aps", reader.aps, "Pa"), ("bps", reader.bps, "1")):
            # 保留原混合坐标系数，使用者可以从 ps 再核算每个层中心的压力。
            variable = output.createVariable(name, "f8", ("level",))
            variable[:] = values
            variable.units = units
        for name, values in (("orography", reader.orography), ("areoid_radius", reader.areoid_radius)):
            # 静态地形仅写一次；排序与动态场一致，防止南北倒置或经度错位。
            variable = output.createVariable(name, "f4", ("lat", "lon"))
            variable[:] = orient_lat_lon(values, reader.latitude, reader.longitude, "native")[0]
            variable.units = "m"
        for index, point in enumerate(points):
            fields, metadata = build_state(reader, point, parameters, variability)
            if variability is not None:
                # EOF 相位会随时刻变化，拆成时间变量；其余固定配置留在全局 JSON 属性中。
                sample_metadata = json.loads(metadata["variability_metadata"])
                sample_phase = sample_metadata.pop("sample_phase_sols", point.elapsed_sols)
                metadata["variability_metadata"] = json.dumps(sample_metadata, sort_keys=True)
                if index == 0:
                    variable = output.createVariable("eof_sample_phase_sols", "f8", ("time",))
                    variable.units = "sol"
                output["eof_sample_phase_sols"][index] = sample_phase
            if index == 0:
                # 首个样本决定输出变量及维度，后续样本复用相同定义而只追加该时刻数值。
                output.setncatts(metadata)
                for name, values in fields.items():
                    dimensions = ("time", "lat", "lon")
                    if values.ndim == 3:
                        # 界面场有 nlevel+1 层，普通三维场有 nlevel 层，二者不能共用层维。
                        vertical = "interface" if name in ("pressure_interface", "height_interface") else "level"
                        dimensions = ("time", vertical, "lat", "lon")
                    # 浮点场以 float32 压缩保存；每个压缩块是一整张单时刻场，便于逐时读写。
                    variable = output.createVariable(name, "f4", dimensions, zlib=True,
                        complevel=1, shuffle=True, chunksizes=(1,) + values.shape, fill_value=np.nan)
                    unit_name = name.removeprefix("increment_")
                    # 增量沿用被扰动物理量的单位；柱质量为 kg/m²，其余未列出的诊断默认无量纲。
                    variable.units = UNITS.get(unit_name, "kg m-2" if name.startswith("column_dust_mass") else "1")
                    if name.endswith("ice_roundoff_correction"):
                        variable.units = "kg m-2"
                        variable.long_name = "positive tiny-negative correction; original=ice-this_correction"
                    if name.startswith("rms_"):
                        # 属性再次强调 RMS 无符号，不能直接从月均值中相减以得到每日状态。
                        variable.long_name = "unsigned RMS variability scale; not a signed realization"
                    if name.startswith("tau_"):
                        variable.long_name = name.replace("_", " ")
                    if name == "eta_diagnostic":
                        variable.diagnostic_only = "true; not native dry hydrostatic WRF eta"
                    if name == "tau_wrf_native700":
                        # 700 Pa 只属于这个历史辐射诊断字段，热力基态 P0 仍然固定为 610 Pa。
                        variable.reference_pressure_pa = 700.0
                        variable.comment = "historical radiation diagnostic, independent of thermodynamic P0=610 Pa"
            for name, values in fields.items():
                # 一次只持有和写入一个时刻；离开外层 with 后关闭并发布完整文件。
                output[name][index] = values
    return path


def export_intermediate(path, reader, point, pressure_levels, variable_table,
    parameters=None, extrapolation="missing", variability=None):
    """Optional WPS v5 export; missing outside columns is rejected, not hidden.

    Endpoint extrapolation is explicit and remains forbidden below terrain.
    A real-domain below-terrain extension needs a separate metgrid policy.
    Native-column products retain source terrain and never use a fill policy.
    """
    # 本函数另写 WPS v5 中间文件；先完成全部场的准备检查，再打开二进制写出器。
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    fields, metadata = build_state(reader, point, parameters, variability)
    # 变量表和 FieldSpec 复用 converter；只导出 enabled 行，不在这里生成缺省尘埃场。
    specs = [field for field in read_field_table(variable_table) if field.enabled]
    if any(spec.aim_v in ("TRC01", "TRC02") for spec in specs):
        # 尘埃需要目标层质量守恒映射及 real 接收验证；普通等压面点插值不能替代这条链。
        message = "dust export requires an audited conservative target-layer/real bridge; point sampling is refused"
        raise ValueError(message)
    levels = np.asarray(pressure_levels, dtype=float)
    # WPS 输出层必须是一维、正值且由下向上压力严格递减；重复层会使后续垂直解释不明确。
    if levels.ndim != 1 or levels.size < 2 or np.any(levels <= 0) or np.any(np.diff(levels) >= 0):
        raise ValueError("WPS pressure levels must be positive, unique and descending")
    lats, lons = np.sort(reader.latitude), np.sort(reader.longitude)
    # 此处沿用 MCD 原生水平分辨率，只提供网格描述；3389.92 是火星半径，单位 km。
    grid = RegularLatLonGrid(lats[0], lats[-1], len(lats), lons[0], lons[-1], len(lons), 3389.92)
    prepared = []
    for spec in specs:
        # 静态地形直接排序；二维场直接取值；压力场显式展开为每个格点相同的目标等压面。
        if spec.src_v == "orography":
            values = orient_lat_lon(reader.orography, reader.latitude, reader.longitude, "native")[0]
        elif spec.kind == "2d":
            values = fields[spec.src_v]
        elif spec.src_v == "__pressure__":
            values = np.broadcast_to(levels[:, None, None], (len(levels), len(lats), len(lons)))
        else:
            # 其余三维字段从源混合层插到目标压力层；缺测/端点外推由显式选项控制。
            # surface_pressure_pa 限制地下层，endpoint 也不能越过地形往地下外推。
            values = interpolate_to_pressure(fields[spec.src_v], fields["pressure"], levels,
                surface_pressure_pa=fields["ps"], extrapolation=extrapolation)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{spec.aim_v}: pressure grid needs explicit extrapolation; no file written")
        prepared.append((spec, values))
    with WrfIntermediateWriter(path, grid, point.hdate, "MCD6.2 background") as writer:
        for spec, values in prepared:
            # 三维场逐压力层输出二维 slab；二维场临时加长度为 1 的维度复用相同写出循环。
            for level, slab in enumerate(values if values.ndim == 3 else values[None]):
                # 三维 xlvl 单位 Pa；二维使用 WPS 的地表层特殊编码 200100。
                xlvl = float(levels[level]) if values.ndim == 3 else 200100.0
                writer.write_field(spec.aim_v, slab, spec.units, spec.desc,
                    xlvl=xlvl, is_wind_earth_rel=int(spec.aim_v in ("UU", "VV")))
    # 重新读回检查日期、水平尺寸和必需字段；这验证文件结构，不证明 WPS/real 已消费字段。
    required = {spec.aim_v for spec in specs if spec.required}
    validate_intermediate_file(path, point.hdate, (len(lats), len(lons)), required)
    return path
