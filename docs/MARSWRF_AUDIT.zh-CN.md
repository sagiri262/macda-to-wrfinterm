# MarsWRF 时间、物理过程与 MACDA 转换审计

本审计基于当前工作区中以下源码和实际数据文件：

- `marswrf/WRFV3/share` 和 `marswrf/WRFV3/external/esmf_time_f90`
- `marswrf/WRFV3/phys`
- `DATA/macda-mro-mcs/mro-mcs-reanalysis_mars_MY28SOY507_MY28SOY537_v2-0.nc`
- `marswrf/em_global_mars/wrfinput_d01`
- `marswrf/em_global_mars/wrfout_d01_0003-00001_00:00:00`
- `marswrf/em_global_mars/wrfrst/wrfrst_d01_0001-00520_00:00:00`

## 1. MarsWRF 如何表示时间

MarsWRF 有两种编译时确定的时间模式。`configure` 中的 `mars` 选择传统模式，`mars24`/`marsda` 选择 `MARS24_TIMING` 模式。

### 已检查文件使用的传统模式

`WRF_PLANET` 把 WRF 日期字符串从公历格式 `YYYY-MM-DD_HH:MM:SS` 改为行星格式 `YYYY-DDDDD_HH:MM:SS`：

- `share/module_date_time.F:203-210, 228-235` 把第 6 到 10 列作为一个五位 day/sol 字段读取，并把月份设为零。
- `share/module_date_time.F:251-264, 541-560` 用 `PLANET_YEAR` 检查 sol 是否有效。
- `share/module_date_time.F:668-722, 742-765` 处理 sol 跨行星年进位，并用 `I5.5` 格式输出。
- `external/esmf_time_f90/Meat.F90:202-209` 把一个火星年设为 669 个模型日。
- `external/esmf_time_f90/ESMF_Time.F90:751-756` 输出同样的五位行星日期字符串。

一个模型日仍用 86400 个时钟单位表示，这样已有的 WRF alarm 和日期运算代码可以继续使用。但这不表示物理上的一个 sol 只有 86400 SI 秒。`frame/module_driver_constants.F:95-101` 定义 `P2SI=1.027491252`，即每个 sol 约为 88775.2 SI 秒；辐射和陆面物理在需要 SI 时间时使用该比例。

`share/set_timekeeping.F` 根据 namelist 中的年、日、时等值建立 WRFU 时钟，并设置 history、restart、boundary 和 auxiliary alarm。`module_bc_time_utilities.F` 只比较这些时钟对象，以决定读取侧边界的时机。`wrf_timeseries.F` 负责输出站点诊断结果，不执行地球时间到火星时间的转换。

### MARS24 模式

启用 `MARS24_TIMING` 后，`PLANET_YEAR` 扩大到 100000，使五位 day 字段能够保存连续的 Mars Solar Date（MSD），而不是每 669 sol 重复一次的年内 sol。`share/module_mars24.F` 实现：

- UTC Julian day 到 TT，并处理闰秒偏移（`:25-70`）
- TT/J2000 偏移（`:72-77`）
- J2000 到太阳黄经 Ls（`:79-150`）
- MSD 与 J2000 之间的转换（`:167-179`）
- MTC、地方平太阳时和地方真太阳时（`:181-225`）

`share/module_planet_utilities.F:4-19` 和 `phys/module_radiation_driver.F:2513-2630` 使用这些函数计算 Ls、太阳赤纬、日心距离和时差。这是源码中唯一进行地球 UTC/J2000 到火星时间天文转换的路径。普通 NetCDF I/O 路径本身不会在火星日期旁另存一份地球日期。

### 实际 WRF 文件证据

已检查文件使用传统模式，而不是 MARS24 模式：

| 文件 | `Times` | 关键属性 |
|---|---|---|
| `wrfinput_d01` | `0001-00001_00:00:00` | `PLANET_YEAR=669`，`P2SI=1.027491` |
| `wrfout_d01_0003-00001_00:00:00` | `0003-00001_00:00:00` | 同上 |
| `wrfrst_d01_0001-00520_00:00:00` | `0001-00520_00:00:00` | 同上 |

`share/output_wrf.F:709-762` 把行星常数写入 NetCDF。`share/input_wrf.F:174-217` 按行星的 year/sol 布局读取 `SIMULATION_START_DATE`。所以 `START_DATE`、`SIMULATION_START_DATE`、`Times` 变量和输出文件名都采用同一种 MarsWRF 字符串格式。

## 2. MACDA 时间对齐

输入文件包含 360 个时次，间隔为两个火星小时。其相互独立的时间坐标为：

- `time=3180.08333333333 .. 3210`，单位为从 MY24 sol 1 的 00:00 MTC 起算的 sol
- `MY_Ls=28`
- `Ls=264.283 .. 283.5132` 度
- `Mars_date=+0028-10-07T02:00:00A .. +0028-10-37T00:00:00A`

文档定义的 MACDA sol 日历按 `669,668,669,668,669` 年长度循环。668-sol 年中各月长度为：

```text
56, 55, 56, 55, 56, 56, 55, 56, 55, 56, 56, 56
```

669-sol 年的第 12 月增加一个 sol。前 9 个月共 500 sol，因此 MY28 第 10 月第 7 sol 是年内第 507 sol。转换结果为：

```text
+0028-10-07T02:00:00A
  -> MY 28，SOY 507，02:00:00 MTC
  -> MarsWRF 0028-00507_02:00:00
  -> HDATE   0028-00507_02:00:00.0000
```

转换器还独立重建连续时间：

```text
MY24 668 + MY25 669 + MY26 669 + MY27 668
+ (SOY507 - 1) + 2/24 = 3180.08333333333 sol
```

如果重建值与 NetCDF `time` 坐标不一致，转换立即报错。`XFCST` 使用连续 `time` 相对本次所选首条记录的差值，再乘以 24 个火星小时，因此在所选强迫序列内保持连续。

旧转换器只删除 `+`、`A` 并把 `T` 换成下划线，会把 `0028-10-07_02:00:00` 交给一个要求第 6 到 10 列为单个 sol 序号的软件。这不是有效的 MarsWRF 日期。

## 3. MACDA 到 WPS/WRF 字段

当前生效的字段映射和转换为：

| MACDA | Intermediate | 处理 | 下游用途 |
|---|---|---|---|
| `temp` | `TT` | sigma 层转压力层 | 最终 `TT` |
| `uwind` | `UU` | sigma 层转压力层 | 最终交错网格 `UU` |
| `vwind` | `VV` | sigma 层转压力层 | 最终交错网格 `VV` |
| 派生 | `PRESSURE` | 每个压力层写常数平面 | metgrid 派生 `PRES` |
| `geop` | `GHT` | 除以 3.72 m/s2 | 最终 `GHT` |
| `psurf` | `PSFC` | 原值，Pa | 地表气压 |
| `tsurf` | `SKINTEMP` | 原值，K | 地表温度 |
| 零场 | `SPECHUMD` | 压力层零场 | 干燥初始状态 |
| 零场 | `QV` | 压力层零场 | 水汽为零的状态 |
| `coldust` | `TAU_OD2D` | `coldust/psurf*700 Pa` | 可选火星沙尘光学厚度 |
| `co2ice` | `CO2ICE` | 原值 | 可选地表 CO2 冰 |

MACDA 气压严格按 `p(k,j,i)=psurf(j,i)*lev(k)` 计算。三维场在对数气压坐标中线性插值到配置的火星压力层。源数据纬度从北到南排列，处理时反转成升序；经度从 `[-180,175]` 归一化到 `[0,355]`；风场标记为相对地理坐标。

`omega` 是气压垂直速度，不能直接复制为 WRF 的几何垂直速度 `W`。`swflux` 和 `lwflux` 是 MarsWRF 会重新计算的地表诊断通量。`dustmmr` 是一个 1.5 微米粒径分布，而当前 MarsWRF 配置使用两个沙尘 tracer/粒径档；将其分配给 `TRC01/TRC02` 必须有明确的科学拆分方案，因此没有擅自映射。

## 4. `phys` 中与火星有关的文件

火星专用实现：

- `module_ra_mars_common.F`：沙尘廓线、MCD/MGS/Viking/MCS/TES 沙尘，以及气溶胶短波和长波加热（`:47-1449`）
- `module_ra_mars_kdm.F`：火星相关 k 分布辐射（`:142-2944`）
- `module_ra_mars_wbm.F`：WBM 可见光/红外辐射（`:15-1005`）
- `module_ra_mars_burk.F`：Burke 短波加热（`:11-377`）
- `module_ra_mars_uv.F`：紫外加热（`:110-494`）
- `module_ra_houben.F`：简化火星牛顿冷却（`:28-145`）
- `module_mp_mars_common.F`：火星粒子沉降和 Stokes 支持（`:18-454`）
- `module_mp_mars_co2_simple.F`：CO2 凝结/升华（`:11-305`）
- `module_mp_mars_h2o_simple.F`：H2O 云微物理（`:11-164`）
- `module_mp_mars_basudustlifting.F`：Basu 扬尘（`:13-234`）
- `module_mp_sedim_dust.F`、`module_mp_sedim_water.F`：沙尘/水沉降
- `module_sf_mars_cendustlifting.F`：双粒径档、单粒径档、可变 tau 和指定沙尘注入（`:22-787`）
- `module_mp_mars_chem.F`：被动化学 tracer 倾向（`:14-408`）
- `module_sf_planet_simple.F`：火星地表/地下能量求解与初始化（`:12-812`，火星选择重点位于 `:416-530`）

包含大量 `WRF_MARS` 分支的共用驱动：

- `module_radiation_driver.F`：选择火星辐射方案，计算轨道/Ls、沙尘和云光学廓线（`:119-209, 616-670, 1434-1541, 1875-2052, 2334-2408, 2513-2946`）
- `module_microphysics_driver.F`：连接 CO2、H2O、沙尘沉降、扬尘和 tracer（`:52-142, 309-386, 466-509, 1149-1543`）
- `module_surface_driver.F`：传递火星冰、沙尘、热惯量字段，并调用行星地面物理（`:126-212, 530-560, 1160-1301, 2148-2220, 3662-3810`）
- `module_pbl_driver.F`：传递火星 tracer 和行星边界层参数（`:77-136, 570-688, 871-878, 1067-1148, 1244-1251, 1591-1598`）
- `module_physics_init.F`：初始化火星辐射、地表、沙尘、CO2/H2O 和 tracer 方案（`:33-38, 253-259, 572-633, 1085-1093, 1289-1502, 1647-1654, 1997-2077`）
- `module_physics_addtendc.F`：把火星标量/tracer 倾向耦合到动力过程（`:40-44, 93-100, 145-149, 237-282, 365-414, 479-528, 695-744`）
- `module_bl_mrf.F`、`module_bl_myjpbl.F`、`module_bl_ysu.F`：CO2/火星热力学和附加 tracer 输送
- `module_sf_sfclay.F` 和 `_3012.F`：火星气体常数与近地层处理

## 5. 已验证范围和下游遗留问题

当前 Python 路径已经能够完成 MACDA 到 WRF intermediate 的转换和自校验。WPS 4.6 的 `rd_intermediate.exe` 也成功读取全部 137 条记录，并报告预期的 72 x 36 火星网格、3389.92 km 半径、数据源 `MACDA` 和行星 HDATE。

但是，已检查的 `marswrf/WPS` 日期工具没有包含 WRFV3 中的 `WRF_PLANET` 五位 sol 分支，并且本地没有可供端到端测试的 `metgrid.exe` 或 `real.exe`。因此，在 `metgrid` 能够按时间调度正确的 Mars HDATE 之前，可能仍需要把 WRF 的行星日期逻辑移植到这份旧 WPS 源码中。样例 namelist 记录的是正确目标日期，但不会掩盖这个源码层面的限制。

可选的 `TAU_OD2D` 和 `CO2ICE` 记录还需要使用所提供的 `METGRID.TBL.MARS.additions` 条目。运行 geogrid/real 仍然需要火星静态地理数据集。

最后，传统 MarsWRF 固定把每一年都视为 669 sol，而 MACDA sol 日历包含 668-sol 年。本次要求处理的 MY28 文件不跨年。对于跨年强迫数据，应统一修改 WRF/WPS 日历，或者使用经过正确 epoch 对齐的 MARS24/MSD 编译；不能悄悄依赖传统模式的固定年长运算。

