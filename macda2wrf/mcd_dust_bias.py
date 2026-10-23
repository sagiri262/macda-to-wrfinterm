"""Fixed-bin dust sensitivities, with optical and settling assumptions separated.

Nothing here changes dust mass, model configuration or a radiation calculation.
Optical stages are additive counterfactual diagnostics in the stated order.
Settling references are specified closures, not observations or weather errors.
"""

# 本模块量化固定双粒径表示带来的敏感性，不修改 MCD 质量场或 MarsWRF 配置。
# spherical_area_weights 给出精确球面面积权重；dust_bias_columns 返回单时刻的二维柱诊断。
# 输出尽量保留可累加的分子量，空间/时间汇总后再求比，避免直接平均局地比值的偏差。

from __future__ import annotations

import numpy as np

# 质量分 bin 和基本光学诊断复用 mcd_dust；速度公式在 mcd_settling 中独立实现。
# 辐射权重、辐射下限只进入光学对照，不参与物理质量或沉降速度的计算。
from .mcd_dust import DustBinParameters, column_dust_diagnostics, pressure_layer_thickness
from .mcd_settling import mcd_deposition_velocity, wrf_lognormal_mass_velocity, wrf_stokes_velocity


def spherical_area_weights(latitude, longitude):
    """Fractional sphere area of a regular global grid, including finite polar caps."""
    # 输入是一维格点中心坐标，输出为 (lat, lon) 面积分数，总和为 1，行星半径因此约去。
    latitude = np.asarray(latitude, dtype=float)
    longitude = np.asarray(longitude, dtype=float)
    # 下面的界面推断只适用于规则全球网格；区域网格、缺测坐标或倒序坐标直接拒绝。
    if latitude.ndim != 1 or longitude.ndim != 1 or min(latitude.size, longitude.size) < 2:
        raise ValueError("latitude and longitude must be one-dimensional global coordinates")
    if not np.all(np.isfinite(latitude)) or not np.all(np.isfinite(longitude)):
        raise ValueError("nonfinite geographic coordinates")
    if np.any(np.diff(latitude) <= 0) or np.any(np.diff(longitude) <= 0):
        raise ValueError("geographic coordinates must be strictly ascending")
    if latitude[0] < -90 or latitude[-1] > 90:
        raise ValueError("invalid latitude")
    # 等间距检查分别作用于纬度与经度，允许的绝对舍入差为 2e-5 度。
    latitude_step = np.diff(latitude)
    longitude_step = np.diff(longitude)
    if not np.allclose(latitude_step, latitude_step[0], rtol=0, atol=2.0e-5):
        raise ValueError("expected regular global latitude grid")
    if not np.allclose(longitude_step, longitude_step[0], rtol=0, atol=2.0e-5):
        raise ValueError("expected regular global longitude grid")
    # Nlon*dlon=360 保证完整周期且不重复计入首尾经线；否则全球面积归一化不成立。
    if not np.isclose(longitude_step[0] * longitude.size, 360, rtol=0, atol=2.0e-5):
        raise ValueError("longitude must cover a full nonduplicated cyclic grid")
    if latitude[0] > -90 + latitude_step[0] / 2 + 2.0e-5:
        raise ValueError("latitude does not reach the south polar cap")
    if latitude[-1] < 90 - latitude_step[0] / 2 - 2.0e-5:
        raise ValueError("latitude does not reach the north polar cap")
    # 纬带边界取相邻中心中点，最外端固定到两极；极点格点也对应有限极帽面积。
    bounds = np.concatenate(([-90.0], (latitude[:-1] + latitude[1:]) / 2, [90.0]))
    # 球面纬带面积/全球面积=(sin(phi_north)-sin(phi_south))/2，不能只取 cos(极点)=0。
    bands = np.diff(np.sin(np.deg2rad(bounds))) / 2
    # 等间距经线把每个纬带平分为 Nlon 份；copy 将广播视图变为独立二维权重数组。
    return np.broadcast_to(bands[:, None] / longitude.size, (latitude.size, longitude.size)).copy()


def dust_bias_columns(dustq, reff_m, pressure_interfaces,
    pressure_center, temperature_k, density_kg_m3, parameters=None):
    """Return 2-D extensive diagnostics from one (level, lat, lon) state.

    All optical stages use 610-Pa normalization. Source optics is a visible
    approximation; WRF Q values do not by themselves prove spectral identity.
    Size-only keeps
    source Q/rho; density, efficiency, radiation weights and floor then change
    in that order. Sequential contributions depend on order, but sum exactly.
    Column mass*velocity is a mass-weighted speed numerator, not a deposition
    flux. Only the bottom-level rho*q*v proxy has kg/m2/s flux units; neither
    proxy reproduces the finite-step Van Leer deposition/transport solver.
    """
    # dustq 是 kg/kg，reff_m 是 m，pressure_* 是 Pa，temperature_k 是 K，density 是 kg/m³。
    # 大气场均为 (level,lat,lon)，界面比层中心多一层；垂直方向固定为地表到模式顶。
    params = parameters or DustBinParameters()
    # 归因实验要求两套重力相同，否则质量积分差异也混入光学阶段，失去受控比较含义。
    if not np.isclose(params.gravity_m_s2, params.wrf_gravity_m_s2, rtol=0, atol=1.0e-12):
        raise ValueError("optical attribution requires identical source and WRF gravity")
    mass = np.asarray(dustq, dtype=float)
    radius = np.asarray(reff_m, dtype=float)
    pressure = np.asarray(pressure_center, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    density = np.asarray(density_kg_m3, dtype=float)
    # 转双精度后核对形状、有限性与物理正值；缺测或非法输入不参与全球平均。
    if mass.ndim != 3 or any(item.shape != mass.shape for item in (radius, pressure, temperature, density)):
        raise ValueError("expected identical (level, lat, lon) input shapes")
    if np.any(~np.isfinite(mass)) or np.any(mass < 0):
        raise ValueError("dust mass must be finite and nonnegative")
    for name, values in (("pressure", pressure), ("temperature", temperature), ("density", density)):
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(f"{name} must be finite and positive")
    if np.any((mass > 0) & (~np.isfinite(radius) | (radius <= 0))):
        raise ValueError("nonzero dust needs finite positive radius")
    interfaces = np.asarray(pressure_interfaces, dtype=float)
    thickness = pressure_layer_thickness(interfaces, mass.shape[0])
    # thickness 是正 dp，不是几何厚度；必须与每层 dustq 同形，并且界面向上压力下降。
    if thickness.shape != mass.shape or np.any(np.diff(interfaces, axis=0) >= 0):
        raise ValueError("require bottom-to-top pressure interfaces matching dust columns")
    if np.any(pressure >= interfaces[:-1]) or np.any(pressure <= interfaces[1:]):
        raise ValueError("pressure centers must lie strictly inside layer interfaces")
    # 底层函数按质量守恒分成两个物理 bin，并记录粒径超范围截取和辐射假设带来的残差。
    dust = column_dust_diagnostics(mass, radius, interfaces, params)
    first, second = dust["q_bin1"], dust["q_bin2"]
    # 零尘处的半径本无定义，诊断运算临时取 r1 防止 0*NaN；其质量权重仍为零。
    represented = np.where(mass > 0, dust["reff_matched"], params.radius1_m)
    radius = np.where(mass > 0, radius, params.radius1_m)
    # 静力柱积分：每层空气柱质量 dp/g，乘 dustq 得 kg dust/m²。
    layer_mass = mass * thickness / params.gravity_m_s2
    # 610/ps 只用于柱光学厚度归一化；这里 interfaces[0] 就是地表压力。
    reference = 610.0 / interfaces[0]
    integral = thickness / params.gravity_m_s2 * reference
    # 球形粒子单位质量消光系数为 3Q/(4*rho_particle*r)，因此核心分子是各 bin 的 q/r。
    first_moment = first / params.radius1_m
    second_moment = second / params.radius2_m
    # 第一对照阶段只把粒径改成固定双 bin，仍保留源 Q 和源颗粒密度。
    size_only = np.sum((first_moment + second_moment) * integral, axis=0)
    size_only *= 0.75 * params.mcd_qext_visible / params.mcd_particle_density
    # 第二阶段只把颗粒密度改为 WRF 辐射密度；默认值对应 2500/2650 的缩放。
    density_only = size_only * params.mcd_particle_density / params.wrf_radiation_density
    # 第三阶段再采用 WRF 的各 bin 效率 Q；尚不乘辐射权重，也不施加下限。
    no_weight = dust["tau_bins_visible610_unweighted"]
    # 第四阶段加入 weight1/weight2，但仍无 floor；不是按权重重新分配或放大物理 dustq。
    weighted = params.qext1 * params.weight1 * first_moment + params.qext2 * params.weight2 * second_moment
    weighted = np.sum(weighted * thickness / params.wrf_gravity_m_s2 * reference, axis=0)
    weighted *= 0.75 / params.wrf_radiation_density
    # 保存二维柱分子：低于 r1/高于 r2 的柱质量用于计算截取质量比例。
    # 半径矩是 sum(layer_mass*radius)，除以总柱尘质量后才得到质量加权平均半径。
    # optical 阶段最终值 tau_native610 额外包含辐射下限；全部阶段统一归一到 610 Pa。
    result = {"column_mass": dust["column_dust_mass"],
        "column_mass_low": np.sum(np.where(radius < params.radius1_m, layer_mass, 0), axis=0),
        "column_mass_high": np.sum(np.where(radius > params.radius2_m, layer_mass, 0), axis=0),
        "radius_requested_mass_moment": np.sum(layer_mass * radius, axis=0),
        "radius_matched_mass_moment": np.sum(layer_mass * represented, axis=0),
        "radius_absolute_error_mass_moment": np.sum(layer_mass * np.abs(represented - radius), axis=0),
        "tau_source_spherical610": dust["tau_visible610_spherical"], "tau_size_only610": size_only,
        "tau_density610": density_only, "tau_efficiency610": no_weight,
        "tau_weight610": weighted, "tau_native610": dust["tau_wrf_native610"]}
    # 依既定顺序相减得到五项增量；望远镜求和严格等于最终值减起始值。
    # 单项归因依赖替换顺序，不表示每个因素在任意顺序下都有同样的独立贡献。
    stages = [result[name] for name in ("tau_source_spherical610", "tau_size_only610", "tau_density610",
        "tau_efficiency610", "tau_weight610", "tau_native610")]
    for index, label in enumerate(("size", "density", "efficiency", "weight", "floor")):
        result[f"tau_delta_{label}610"] = stages[index + 1] - stages[index]
    result["tau_size_absolute_error610"] = np.abs(stages[1] - stages[0])
    # 绝对差和平方差独立保留，汇总后计算 MAE/RMSE；它们是方案差异，非观测真值误差。
    result["tau_size_squared_error610"] = (stages[1] - stages[0])**2
    # 沉降对照统一采用 namelist 颗粒密度和重力，不使用光学计算的 2650 辐射密度。
    arguments = {"particle_density": params.wrf_namelist_density, "gravity_m_s2": params.wrf_gravity_m_s2}
    # 先在相同层中心 p/T 上算 r1、r2 的向下终端速度，单独考察粒径表示造成的差别。
    velocity1 = wrf_stokes_velocity(pressure, temperature, params.radius1_m, **arguments)
    velocity2 = wrf_stokes_velocity(pressure, temperature, params.radius2_m, **arguments)
    bin_speed_mass = (first * velocity1 + second * velocity2) * thickness / params.gravity_m_s2
    result["speed_bins_center_mass_moment"] = np.sum(bin_speed_mass, axis=0)
    # 再改用各层下界面压力并保留层温，贴近 WRF 速度函数的调用约定；速度变量在此被覆盖。
    velocity1 = wrf_stokes_velocity(interfaces[:-1], temperature, params.radius1_m, **arguments)
    velocity2 = wrf_stokes_velocity(interfaces[:-1], temperature, params.radius2_m, **arguments)
    native_speed_mass = (first * velocity1 + second * velocity2) * thickness / params.gravity_m_s2
    result["speed_bins_interface_mass_moment"] = np.sum(native_speed_mass, axis=0)
    # 单一 reff 与假设 veff=0.5 对数正态分布都在层中心 p/T 上比较，不能当成观测速度。
    # 有效半径本身不足以确定粒径分布，宽分布速度是额外闭合假设的结果。
    mono = wrf_stokes_velocity(pressure, temperature, radius, **arguments)
    broad = wrf_lognormal_mass_velocity(pressure, temperature, radius, **arguments)
    result["speed_monodisperse_mass_moment"] = np.sum(layer_mass * mono, axis=0)
    result["speed_lognormal_v05_mass_moment"] = np.sum(layer_mass * broad, axis=0)
    # 上述 speed_*_mass_moment 单位 kg/(m*s)，是平均速度分子，不是穿过地面的通量。
    # 此处底层 rho*(q1*v1+q2*v2) 才有 kg/(m²*s) 通量单位，且用的是下界面 bin 速度。
    result["surface_flux_bins_proxy"] = density[0] * (first[0] * velocity1[0] + second[0] * velocity2[0])
    mcd_speed = mcd_deposition_velocity(interfaces[0], temperature[0], radius[0], params.gravity_m_s2)
    # 两个地表代理量共用 MCD 最低层 rho、dustq，只更换速度配方；不是原生模式沉降输出。
    # 没有做有限时间步的沉降输送、Van Leer 重构或地表质量收支积分。
    result["surface_flux_mcd_speed_recipe_same_q"] = density[0] * mass[0] * mcd_speed
    return result
