"""Stream fixed 0.5/1.2-micron bin diagnostics without rewriting MY36 backgrounds.

Run from the project root: python -m utils.check_mcd_dust_bins. Equal 2h samples
are weighted by exact spherical cell area; radius and speed means additionally
use layer dust mass q*dp/g. Exclude the four wrapped records from seasonal and
annual summaries, but retain their time diagnostics and flags. Never overwrite.
"""

# 诊断入口：读取已保存的未扰动 MY36 背景，重新计算固定双 bin 的光学与沉降敏感性。
# 先核验全年时间覆盖，再逐 sol 文件、逐时刻处理，最后写入一个独立诊断 NetCDF。
# 它不重写背景文件，也不修改 WRF 粒径/namelist；建议从项目根目录以 -m 模式运行。
# 全部统计保留面积、质量和时间权重的定义，避免局地比值平均或季节样本数不同导致误读。

from __future__ import annotations

# argparse 提供命令行入口；参数 dataclass 转 JSON 后与背景元数据及输出一同保留。
import argparse
from dataclasses import asdict
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np

# 分 bin、柱诊断和面积权重直接复用物理模块；new_dataset 保证完整关闭后发布且不覆盖。
from macda2wrf.mcd_dust import DustBinParameters, split_dust_bins
from macda2wrf.mcd_dust_bias import dust_bias_columns, spherical_area_weights
from macda2wrf.mcd_pipeline import new_dataset
from macda2wrf.mcd_reader import read_array
from .check_mcd_reconstruction import check_coverage


# 默认输入是逐 sol 背景目录，输出为另外一个诊断文件；调用者可通过命令行指定其他路径。
DEFAULT_INPUT = "/home/zy/WRF/DATA/mcd_reanalysis/background"
DEFAULT_OUTPUT = "/home/zy/WRF/DATA/mcd_reanalysis/MY36_fixed_dust_bin_diagnostics.nc"
# 每个输入变量都逐字检查单位；tau_visible610_mcd 是源可见光柱厚，不是 9.3 μm 吸收厚度。
INPUT_UNITS = {"dustq": "kg kg-1", "reffdust": "m", "pressure_interface": "Pa",
    "pressure": "Pa", "temperature": "K", "rho": "kg m-3", "tau_visible610_mcd": "1",
    "tau_wrf_native610": "1", "q_bin1": "kg kg-1", "q_bin2": "kg kg-1"}
# 派生比例/平均值的单位另列；原始柱分子量的单位由 field_units 推导。
# 名称中的 bias、MAE、RMSE 表示固定 bin 对照差异，不表示对独立观测真值的验证误差。
METRIC_UNITS = {"clipped_mass_fraction": "1", "low_mass_fraction": "1", "high_mass_fraction": "1",
    "radius_requested_mass_mean": "m", "radius_matched_mass_mean": "m", "radius_bias_mass_mean": "m",
    "radius_mae_mass_mean": "m", "speed_bins_center_mass_mean": "m s-1",
    "speed_bins_interface_mass_mean": "m s-1", "speed_monodisperse_mass_mean": "m s-1",
    "speed_lognormal_v05_mass_mean": "m s-1", "tau_size_relative_change": "1",
    "tau_native_relative_change": "1", "tau_spherical_to_mcd_ratio": "1",
    "speed_bins_to_mono_ratio": "1", "speed_bins_to_lognormal_v05_ratio": "1",
    "surface_flux_bins_to_mcd_ratio": "1", "tau_size_mae": "1", "tau_size_rmse": "1"}


def field_units(name):
    # tau 无量纲；半径分子为柱尘质量*半径，所以单位 kg/m，尚不是 m。
    if name.startswith("tau_"):
        return "1"
    if name.startswith("radius_"):
        return "kg m-1"
    # 速度分子为柱尘质量*速度；只有 surface_flux 是 rho*q*v 的地表通量代理量。
    if name.startswith("speed_"):
        return "kg m-1 s-1"
    if name.startswith("surface_flux_"):
        return "kg m-2 s-1"
    # 其余原始诊断均按柱尘质量处理，单位 kg/m²。
    return "kg m-2"


def ratio(numerator, denominator):
    # 分母没有正物理权重时返回 NaN，避免零尘或零光学厚度被伪装成有效比例。
    return float(numerator / denominator) if denominator > 0 else float("nan")


def summary_metrics(means):
    """Ratios of extensive weighted sums, never averages of local ratios."""
    # means 已是同一范围内的面积/时间加权柱分子；先累加分子分母，再求比。
    # 例如 sum(area*尘质量*半径)/sum(area*尘质量)，不是各格点半径或局地比值的算术平均。
    mass = means["column_mass"]
    requested = ratio(means["radius_requested_mass_moment"], mass)
    matched = ratio(means["radius_matched_mass_moment"], mass)
    source_tau = means["tau_source_spherical610"]
    # 质量超范围比例分成低端、高端及合计；requested/matched 比较原半径与双 bin 能表示的半径。
    result = {"low_mass_fraction": ratio(means["column_mass_low"], mass),
        "high_mass_fraction": ratio(means["column_mass_high"], mass),
        "clipped_mass_fraction": ratio(means["column_mass_low"] + means["column_mass_high"], mass),
        "radius_requested_mass_mean": requested, "radius_matched_mass_mean": matched,
        "radius_bias_mass_mean": matched - requested,
        "radius_mae_mass_mean": ratio(means["radius_absolute_error_mass_moment"], mass),
        # 相对光学变化均以同一源球形近似为分母，不以 WRF 最终值或局地比值为分母。
         "tau_size_relative_change": ratio(means["tau_size_only610"] - source_tau, source_tau),
        "tau_native_relative_change": ratio(means["tau_native610"] - source_tau, source_tau),
        "tau_size_mae": means["tau_size_absolute_error610"],
        "tau_size_rmse": float(np.sqrt(means["tau_size_squared_error610"])),
        # 球形常 Q 近似与源 MCD tau 的比值用于暴露光学假设差异，不强制二者相等。
        "tau_spherical_to_mcd_ratio": ratio(source_tau, means["tau_mcd_source610"]),
        "speed_bins_to_mono_ratio": ratio(means["speed_bins_center_mass_moment"],
            means["speed_monodisperse_mass_moment"]),
        "speed_bins_to_lognormal_v05_ratio": ratio(means["speed_bins_center_mass_moment"],
            means["speed_lognormal_v05_mass_moment"]),
        "surface_flux_bins_to_mcd_ratio": ratio(means["surface_flux_bins_proxy"],
            means["surface_flux_mcd_speed_recipe_same_q"])}
    # 四种速度分子都除以相同柱尘质量，得到可比较的质量加权速度 m/s。
    for label in ("bins_center", "bins_interface", "monodisperse", "lognormal_v05"):
        result[f"speed_{label}_mass_mean"] = ratio(means[f"speed_{label}_mass_moment"], mass)
    return result


def validate_columns(columns, saved_native):
    """Independent invariants plus closure against the already archived native tau."""
    # 五段增量应望远镜闭合到最终 tau 减源球形 tau；双精度重算使用较严格容差。
    reference = columns["tau_source_spherical610"]
    increments = sum(columns[f"tau_delta_{label}610"] for label in ("size", "density", "efficiency", "weight", "floor"))
    np.testing.assert_allclose(increments, columns["tau_native610"] - reference, rtol=1.0e-11, atol=1.0e-12)
    # 本诊断固定默认密度，单独替换密度阶段必须等于前一阶段乘 2500/2650。
    np.testing.assert_allclose(columns["tau_density610"], columns["tau_size_only610"] * 2500 / 2650,
        rtol=1.0e-11, atol=1.0e-12)
    if np.any(columns["tau_delta_floor610"] < -1.0e-12):
        raise ValueError("negative radiation floor increment")
    # 已存背景为 float32，重算为 float64，允许输出量化误差；不允许显著改变已存光学结果。
    np.testing.assert_allclose(columns["tau_native610"], saved_native, rtol=3.0e-6, atol=2.0e-7)


def validate_saved_bins(mass, radius, first, second, parameters):
    """Reject invalid bins; compare splits relative to total mass for float32 inputs.

    Radius rounding near a bin edge can change a tiny bin by a large fraction
    of itself. Use total dust as the scale instead, without accepting NaN,
    negative mass, or a material swap between bins.
    """
    # 按当前物理分 bin 函数重算期望质量，不使用辐射权重去修正物理示踪剂。
    expected = split_dust_bins(mass, radius, parameters)
    for name, saved in (("q_bin1", first), ("q_bin2", second)):
        # 先查形状、有限性、非负性，不能让比较容差接受 NaN 或负物理质量。
        if saved.shape != mass.shape or np.any(~np.isfinite(saved)) or np.any(saved < 0):
            raise ValueError(f"{name}: archived bin must be finite nonnegative physical mass")
        # 误差按总 dustq 归一：粒径 float32 舍入会令近零小 bin 的自身相对误差很大。
        # 用总量作尺度仍能拒绝 bin 间的实质质量交换；tiny 只防止零尘处除零。
        split_relative = np.abs(saved - expected[name]) / np.maximum(mass, np.finfo(float).tiny)
        if np.any(split_relative > 3.0e-6):
            raise ValueError(f"{name}: archived bin differs from mass-conserving radius split")
    # 同时检查两个 bin 相加是否回到原 dustq，而不只检查总量或只检查单个 bin。
    relative = np.abs(first + second - mass) / np.maximum(mass, np.finfo(float).tiny)
    if np.any(relative > 3.0e-6):
        raise ValueError("saved two-bin physical mass does not close")
    return float(np.max(relative))


def write_variable(output, name, values, dimensions, units):
    # 诊断保留 float64；带维度变量启用压缩，全年标量不启用压缩；NaN 继续表示缺测/无定义。
    variable = output.createVariable(name, "f8", dimensions, zlib=bool(dimensions), fill_value=np.nan)
    variable[:] = values
    variable.units = units


def diagnose(input_dir=DEFAULT_INPUT, output_file=DEFAULT_OUTPUT):
    """Validate complete MY36 coverage, stream diagnostics, publish one compact NetCDF."""
    # 先解析绝对路径并拒绝覆盖；只有全部读取、闭合检查和统计成功才进入写出阶段。
    destination = Path(output_file).expanduser().resolve()
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite {destination}")
    # 共用检查器核对完整 MY36 的 8028 个唯一 2h 时刻、日期、Ls 与轨道回绕标志。
    # manifest 仅保存路径和时间元数据，并按文件首时刻排序，不把全年大气场载入内存。
    manifest = check_coverage(input_dir)
    # time_series 存逐时全球数值；seasonal_sums 存各 30° Ls 箱的二维柱分子和。
    time_series, seasonal_sums = {}, {}
    counts = np.zeros(12, dtype=np.int64)
    # longitudes 在这个函数中保存太阳黄经 Ls 时间序列；地理经度另存于 longitude。
    time_values, dates, longitudes, wrapped_flags = [], [], [], []
    latitude, longitude, parameters, area = None, None, None, None
    maximum_native_residual, maximum_mass_relative_residual = 0.0, 0.0
    for file_index, item in enumerate(manifest):
        with Dataset(item["path"]) as dataset:
            # 只接受原未扰动背景，避免把 EOF 扰动候选混入固定 bin 的受控敏感性统计。
            if getattr(dataset, "product_type", "") != "mcd_interpolated_background":
                raise ValueError("diagnostic requires unchanged unperturbed MCD backgrounds")
            current = DustBinParameters(**json.loads(dataset.dust_parameters))
            # 固定到当前默认 0.5/1.2 μm 半径及全部光学默认值，防止异配置产品混用。
            if current != DustBinParameters():
                raise ValueError("this diagnostic is pinned to the current 0.5/1.2um candidate and optical defaults")
            if float(dataset.reference_pressure_pa) != 610.0:
                raise ValueError("expected 610 Pa reference pressure")
            if file_index == 0:
                # 首文件建立网格和面积权重基准；后续文件必须在同一地理格点上。
                parameters = current
                latitude, longitude = read_array(dataset["lat"]), read_array(dataset["lon"])
                area = spherical_area_weights(latitude, longitude)
            for name, expected in (("lat", latitude), ("lon", longitude)):
                np.testing.assert_allclose(read_array(dataset[name]), expected, rtol=0, atol=2.0e-5)
            for name, units in INPUT_UNITS.items():
                if str(getattr(dataset[name], "units", "")) != units:
                    raise ValueError(f"{name}: inconsistent units")
            # 一次读取一个 sol 文件的必要字段；文件关闭后释放句柄，避免全年累积打开文件。
            values = {name: read_array(dataset[name]) for name in INPUT_UNITS}
            for index in range(len(item["time"])):
                # 逐时重算二维柱诊断，底层负责压力层和物理输入的形状/有效性核验。
                columns = dust_bias_columns(values["dustq"][index], values["reffdust"][index],
                    values["pressure_interface"][index], values["pressure"][index],
                    values["temperature"][index], values["rho"][index], parameters)
                columns["tau_mcd_source610"] = values["tau_visible610_mcd"][index]
                # 每个柱诊断必须与面积权重同形且全有限，否则拒绝继续做全球平均。
                for name, column in columns.items():
                    if column.shape != area.shape or not np.all(np.isfinite(column)):
                        raise ValueError(f"{name}: missing or inconsistent diagnostic data")
                validate_columns(columns, values["tau_wrf_native610"][index])
                # 保存全年最坏的光学绝对残差和双 bin 质量相对残差，以便快速检查闭合程度。
                native_residual = np.max(np.abs(columns["tau_native610"] - values["tau_wrf_native610"][index]))
                maximum_native_residual = max(maximum_native_residual, float(native_residual))
                mass_residual = validate_saved_bins(values["dustq"][index], values["reffdust"][index],
                    values["q_bin1"][index], values["q_bin2"][index], parameters)
                maximum_mass_relative_residual = max(maximum_mass_relative_residual, mass_residual)
                # 面积分数总和为 1，所以 sum(column*area) 是全球面积均值，不是乘行星面积的总量。
                # 半径和速度仍保留质量矩，再交给 summary_metrics 除以全球柱尘质量。
                means = {name: float(np.sum(column * area)) for name, column in columns.items()}
                means.update(summary_metrics(means))
                for name, value in means.items():
                    time_series.setdefault(name, []).append(value)
                if not item["wrapped"][index]:
                    # 同一时刻的逐时诊断始终保留；只在季节/全年统计中排除轨道末尾回绕样本。
                    # month 表示 floor(Ls/30°) 的季节箱，不是等时长月份，也不使用地球季节。
                    month = min(int(item["ls"][index] // 30), 11)
                    counts[month] += 1
                    for name, column in columns.items():
                        if name not in seasonal_sums:
                            # 每个原始柱分子只分配 (12,lat,lon) 累加器，不存全年三维大气数组。
                            seasonal_sums[name] = np.zeros((12,) + area.shape)
                        seasonal_sums[name][month] += column
            # 这些时间数组包含回绕记录，以便最终的 global_* 与 wrf_date 一一对应。
            time_values.extend(item["time"])
            dates.extend(item["dates"])
            longitudes.extend(item["ls"])
            wrapped_flags.extend(item["wrapped"])
        if (file_index + 1) % 100 == 0 or file_index + 1 == len(manifest):
            # 进度只发到终端，并立即刷新；本脚本没有主动创建日志文件。
            print(f"checked {file_index + 1}/{len(manifest)} sol files", flush=True)
    if np.any(counts == 0):
        raise ValueError("missing seasonal samples")
    # 每箱除以自己的样本数得季节平均图；全年除以所有有效样本数，不能把 12 箱等权平均。
    # 所有保留时刻都是等间隔 2 个火星小时，因此可用样本数代表相同时间权重。
    monthly = {name: values / counts[:, None, None] for name, values in seasonal_sums.items()}
    annual = {name: np.sum(values, axis=0) / np.sum(counts) for name, values in seasonal_sums.items()}
    annual_global = {name: float(np.sum(values * area)) for name, values in annual.items()}
    # 年均/季均先汇总原始柱分子再求比例，避免把各时刻比例平均当作全年总体比值。
    annual_global.update(summary_metrics(annual_global))
    seasonal_global = {name: np.sum(values * area, axis=(1, 2)) for name, values in monthly.items()}
    for month in range(12):
        # METRIC_UNITS 里的派生量不再次送入分子统计；只用当前季节的原始面积均值求指标。
        metrics = summary_metrics({name: values[month] for name, values in seasonal_global.items()
            if name not in METRIC_UNITS})
        for name, value in metrics.items():
            seasonal_global.setdefault(name, np.full(12, np.nan))[month] = value
    with new_dataset(destination) as output:
        # 写清所有科学限定：非同化产品、610 Pa 归一、常 Q 球形近似、顺序依赖的光学归因。
        output.product_type = "mcd_fixed_dust_bin_sensitivity_diagnostic"
        output.is_reanalysis = "false; no observations assimilated or model integration performed"
        output.reference_pressure_pa = 610.0
        output.background_directory = str(Path(input_dir).resolve())
        output.dust_parameters = json.dumps(asdict(parameters), sort_keys=True)
        output.included_records = int(np.sum(counts))
        output.excluded_wrapped_records = int(np.sum(wrapped_flags))
        output.area_weighting = "exact sine differences of midpoint latitude bounds, finite polar caps"
        output.mass_weighting = "cell area * q * abs(dp) / g; all retained 2h samples equal time weight"
        output.optical_attribution_order = "size; particle density; Qeff; radiation weights; floor"
        # Q 数值不保证波段相同；单分散/对数正态速度是指定闭合，不能声称为观测真值。
        output.optical_baseline = "constant Q2.4/rho2500 spherical visible approximation, not exact MCD radiation"
        output.spectral_warning = "WRF Qeff wavelength equivalence is assumed, not independently verified"
        output.settling_comparators = "WRF single-radius or assumed lognormal veff0.5; not observed true velocities"
        output.settling_viscosity_pa_s = 1.0e-5
        output.mcd_surface_closure = "MCD speed recipe only: rsed2.25reff, mean-free-path shape factor0.5"
        # 通量只对比速度配方且共用 MCD rho/q；重力也统一以控制实验，不冒充原生 MCD 沉降。
        output.surface_flux_warning = "same MCD rho0 and dustq; not original MCD deposition or Van Leer output"
        output.mcd_gravity_policy = "controlled comparison uses WRF3.727, not native MCD3.7257964"
        output.vertical_coordinate = "unchanged native hybrid layers; inferred interfaces; no WRF regridding"
        output.no_model_changes = "WRF source, namelist and archived background files unchanged"
        output.maximum_saved_native_absolute_residual = maximum_native_residual
        output.maximum_saved_bin_mass_relative_residual = maximum_mass_relative_residual
        # 建立逐时、12 季节箱和水平网格四个维；不输出已经柱积分掉的 level 维。
        for name, length in (("time", len(time_values)), ("month", 12),
            ("lat", latitude.size), ("lon", longitude.size)):
            output.createDimension(name, length)
        for name, values, units in (("lat", latitude, "degrees_north"), ("lon", longitude, "degrees_east"),
            ("time", time_values, "s"), ("ls", longitudes, "degree"), ("orbit_wrapped", wrapped_flags, "1"),
            ("month", np.arange(1, 13), "1"), ("month_sample_count", counts, "1")):
            dimension = name
            # ls 和回绕标志依 time；month 与各箱样本数依 month，lat/lon 各依自己的坐标维。
            if name in ("time", "ls", "orbit_wrapped"):
                dimension = "time"
            elif name.startswith("month"):
                dimension = "month"
            write_variable(output, name, values, (dimension,), units)
        # 日期字符串与 time 同长，cell_area_fraction 保存实际使用的二维面积权重以便复核。
        output.createVariable("wrf_date", str, ("time",))[:] = np.asarray(dates, dtype=object)
        write_variable(output, "cell_area_fraction", area, ("lat", "lon"), "1")
        for name, values in time_series.items():
            units = METRIC_UNITS.get(name, field_units(name))
            # 每个指标给出逐时全球序列、12 箱全球序列和全年全球标量，单位在三者中一致。
            write_variable(output, f"global_{name}", values, ("time",), units)
            write_variable(output, f"seasonal_global_{name}", seasonal_global[name], ("month",), units)
            write_variable(output, f"annual_global_{name}", annual_global[name], (), units)
        for name, values in monthly.items():
            # 空间图只保存原始柱分子量的季节均值和年均值，派生比值不在格点上先行混平均。
            write_variable(output, f"monthly_mean_{name}", values, ("month", "lat", "lon"), field_units(name))
            write_variable(output, f"annual_mean_{name}", annual[name], ("lat", "lon"), field_units(name))
    # 完整文件已关闭发布后才打印结果位置、纳入/排除样本数和年均指标；不另写报告文件。
    print(json.dumps({"output": str(destination), "included_records": int(np.sum(counts)),
        "excluded_wrapped_records": int(np.sum(wrapped_flags)), "annual": annual_global}, sort_keys=True), flush=True)
    return destination


def main():
    # 模块命令行只接收输入目录和独立输出路径；计算细节固定于上述受控诊断流程。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default=DEFAULT_INPUT)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    diagnose(arguments.input_dir, arguments.output)


if __name__ == "__main__":
    # import 本模块只提供函数；直接执行或 python -m 时才读取背景并写诊断产品。
    main()
