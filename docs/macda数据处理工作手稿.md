1、macda的数据一共有19个变量，所有关于该再分析数据的文档子和资料都被存放在 ./DOCS 目录下
要了解的东西包括：文件的格式（HDF5），变量名、物理意义、变量的值，

2、对于cmip6-to-wrfinterm目录，conf目录下初始化数据，db目录下存放Vtable，要把MACDA的变量名称对应到WRF的输入变量


**/home/zy/WRF/marswrf/WRFV3/test/em_global_mars** 的 namelist.input、运行脚本、样例 wrfinput 头信息和嵌套案例文件

**/home/zy/WRF/marswrf/WRFV3/Registry** 的 Registry.EM_PLANET、registry.dimspec 等核心变量定义
module_initialize_global_mars.F、module_initialize_real.F、input_wrf.F、module_optional_input.F 中和 ideal/real/WPS 输入相关的逻辑。

MACDA 文件 **/home/zy/WRF/DATA/macda-mro-mcs/mro-mcs-reanalysis_mars_MY28SOY507_MY28SOY537_v2-0.nc** 的 NetCDF 变量、单位、维度和 sigma 垂直坐标说明。

后续任务的核心：要把 MACDA 的 **temp/uwind/vwind/geop/psurf/tsurf/co2ice/coldust/dustmmr/swflux/lwflux** 等变量，转换成 MarsWRF/WPS/real.exe 能接受的输入链路；这会涉及变量名映射、单位换算、sigma 到压力/高度层处理、Mars 静态地表变量补齐、metgrid 可读 intermediate 文件生成，以及 MarsWRF real.exe/wrf.exe 对输入字段和 flags 的要求。


另外先记录两个重要发现：
1、当前 em_global_mars 顶层没有 ideal.exe/real.exe/wrf.exe；只在子目录找到部分 ideal.exe 和 wrf.exe，没有找到 real.exe。
2、现有 Mars case 更像是 ideal.exe -> wrfinput -> wrf.exe 的理想化路径；要走 MACDA -> intermediate -> metgrid -> real.exe -> wrf.exe，需要额外确认或补齐 MarsWRF 的 real-data 编译和输入字段适配。

# 静态数据
有公开数据，有公开方案，但没有可以直接下载后放进 geog_data_path 就能跑 MarsWRF 的完整静态地理包。所以我们要实现geogrid静态数据的路线是：用公开 Mars 数据源自己制作 WPS_GEOG_MARS，再改 GEOGRID.TBL，必要时还要改 real.exe 读取逻辑。
当前 MarsWRF 的 Registry 里确实有 HGT_M、ALBEDO/ALBBCK、EMISS/EMBCK、THC/THCBCK 这些 Mars 静态/地表物理量；但 **/home/zy/WRF/marswrf/WPS/geogrid/GEOGRID.TBL.ARW:4** 仍主要是地球 WPS 默认静态场，例如 HGT_M、LANDUSEF、SOILCTOP、ALBEDO12M 等。所以不能直接套地球 WPS 数据。
公开数据源
## 地形 HGT_M
用 NASA/PDS 的 MOLA MEGDR。这是 Mars Global Surveyor 的 MOLA 全球地形数据，PDS 明确提供 4、16、32、64、128 pixels/degree 分辨率，且包含 topography、radius、areoid、counts 等产品。
链接：PDS MGS MOLA MEGDRs

## 热惯量 THC / THCBCK
用 MGS TES Thermal Inertia Maps。PDS 页面说明这是由 TES thermal bolometer 温度观测反演出来的热惯量图，覆盖 1999-2004，数值范围约 5-5000 J/(m^2 K sqrt(s))。
链接：PDS TES Thermal Inertia Maps
注意：本地 MarsWRF Registry 写的是 THC 单位 J m-1 K-1 s-0.5，而 PDS 热惯量常用物理单位是 J m-2 K-1 s-1/2。这个单位差异后面必须查代码确认，不能盲转。

## 反照率 ALBEDO / ALBBCK
可查 TES bolometric albedo / Mars surface albedo、Mars Express OMEGA 的 bolometric solar hemispherical albedo。PDS TES 标准产品包含 bolometric radiance 和派生 surface properties。
链接：PDS MGS TES
论文关键词可看：Vincendon Mars Express measurements of surface albedo changes 2004 2010
链接：Vincendon et al. arXiv

## 发射率 EMISS / EMBCK
可从 TES emissivity / thermal infrared products 查，但第一版实验也可以先用常数场，例如 0.95 或 1.0，再逐步替换成真实数据。你当前 MarsWRF 理想化初始化里也有手动设定 emiss_0 的逻辑。

## 尘埃光学厚度 / 气溶胶场，后续用
这不一定是 geogrid.exe 的静态地理输入，但 MarsWRF 有 TAU_OD2D 等相关变量，未来会用到。可查 Montabone 的 Mars dust climatology。
链接：Montabone et al. arXiv

## 公开方案
核心方案就是 WPS 官方支持的“自定义 geogrid 静态数据”流程。官方文档说明 geogrid.exe 不只处理地球默认数据，也可以通过 GEOGRID.TBL 插入新的 continuous/categorical 静态场；静态数据要写成 WPS geogrid binary tile，加一个 index 元数据文件。


**实际要做的流程大概是**：
1、下载 MOLA/TES/OMEGA 等公开栅格数据。
2、转成统一 Mars 经纬度网格，通常是 regular_ll，经度 0-360 或 -180-180 要统一。
3、写成 WPS geogrid binary tile，例如 00001-xxxxx.00001-yyyyy。
4、为每个数据目录写 index：projection=regular_ll、type=continuous、dx/dy、known_lat/known_lon、wordsize、scale_factor、units 等。
5、改 GEOGRID.TBL，新增或替换 HGT_M、ALBEDO、THC、EMISS 的数据源路径。
6、跑 geogrid.exe，检查 geo_em.d01.nc 是否真的包含这些字段。
7、再确认 metgrid.exe -> real.exe -> wrfinput_d01 是否保留这些字段；如果 real.exe 不读 THC/EMISS/ALBEDO，就需要改 Registry 或初始化代码。

## 建议相关研究和文档的关键词
WPS geogrid binary format index file
WRF GEOGRID.TBL custom static data
WPS static geographical data regular_ll
MOLA MEGDR Mars topography PDS
Mars Global Surveyor MOLA MEGDR 64 pixels per degree
MGS TES Thermal Inertia Maps PDS
Putzig Mellon Global Thermal Inertia and Surface Properties of Mars
TES bolometric albedo Mars PDS
Mars surface bolometric albedo OMEGA Vincendon
Mars Climate Database surface albedo thermal inertia topography
MarsWRF MOLA thermal inertia albedo emissivity
PlanetWRF Mars surface thermal inertia albedo MOLA

**第一阶段不要追求一次性把所有 Mars 静态场做全。先做最小可验证集：HGT_M、ALBEDO、THC、EMISS。其中 HGT_M 用 MOLA，THC 用 TES Thermal Inertia Maps，ALBEDO/EMISS 可以先用公开图或常数场。然后重点验证 MarsWRF 的 real.exe 是否能把这些字段带进 wrfinput_d01。**

# 最小输入变量
**最小大气变量**
· UU：东西向风，m s-1
· VV：南北向风，m s-1
· TT：温度，K
· PRES：压力，Pa
· GHT：位势高度，m
· RH 或 QV 或 SPECHUMD：湿度相关量，Mars 可先置零或按 MACDA 水汽处理
· PSFC：地表压力，Pa，Mars 这里非常重要
· SKINTEMP：地表皮肤温度，K
· TAVGSFC：日平均近地面温度，K，没有 SKINTEMP 时可兜底
· PMSL：海平面气压，地球概念；Mars 不建议依赖它，最好直接给 PSFC
**METGRID.TBL.ARW** 里 PRES、TT、UU、VV 是 mandatory；real.exe 里还明确检查 TSK 或 TAVGSFC，否则会 fatal。
**静态地理变量**
· XLAT_M、XLONG_M
· HGT_M：地形高度，m
· LANDMASK
· SOILTEMP
· ALBEDO12M
· GREENFRAC
· SCT_DOM、SCB_DOM

但 Mars 需要额外特别注意这些：
· ALBEDO / ALBBCK：反照率
· EMISS / EMBCK：地表发射率
· THC / THCBCK：热惯量，J m-1 K-1 s-0.5
· H2OICE：地表水冰，kg m-2
· CO2ICE：地表 CO2 冰，kg m-2
· TAU_OD2D：归一到 7 mb 的尘埃光学厚度
· TAU_OD：三维光学厚度
· GRD_ICE_PC、GRD_ICE_DP：地下冰比例和深度

问题是：这些 Mars 专用场虽然在 Registry 里有输入/输出标记，但不在 realonly 包里。也就是说，直接把它们写进 met_em 不一定会被 real.exe 正确带入 wrfinput_d01。后面很可能需要改 Registry.EM_PLANET、module_initialize_real.F，或者在 real.exe 后对 wrfinput_d01 做注入。

## wrfinput / wrfout 里的核心变量
**wrfinput_d01 当前维度是**
* `west_east = 72`
* `south_north = 36`
* `bottom_top = 52`
* `bottom_top_stag = 53`
* `soil_layers_stag = 15`
* `west_east_stag = 73`
* `south_north_stag = 37`

**WRF 内部动力变量主要是**
* `U`：X stagger 风，`m s-1`
* `V`：Y stagger 风，`m s-1`
* `W`：Z stagger 垂直风，`m s-1`
* `T`：扰动位温，`K`
* `P`：扰动压力，`Pa`
* `PB`：基态压力，`Pa`
* `PH`：扰动位势，`m2 s-2`
* `PHB`：基态位势，`m2 s-2`
* `MU`：柱积分扰动干空气质量，`Pa`
* `MUB`：柱积分基态干空气质量，`Pa`
* `P_HYD`：静力压力，`Pa`
* `PSFC`：地表气压，`Pa`
* `TSK`：地表皮肤温度，`K`
* `TSLB`：土壤温度，`K`

**Mars 物理相关输出包括**
* `ALBEDO`
* `ALBBCK`
* `EMBCK`
* `THCBCK`
* `EMISS`
* `THC`
* `H2OICE`
* `CO2ICE`
* `TAU_OD2D`
* `TAU_CL2D`
* `TAU_OD`
* `TAU_CL`
* `QV_COLUMN`
* `QI_COLUMN`
* `TRC01`、`TRC02`：两个尘埃 bin 混合比，`kg kg-1`

## 最需要注意的地方

1. **不要把 MACDA 变量直接对应到 `wrfinput` 内部变量。**
   MACDA 到 WPS intermediate 应该写 `UU/VV/TT/PRES/GHT/PSFC/SKINTEMP` 等，`U/V/T/P/PB/PHB` 是 `real.exe` 之后才生成的内部变量。
2. **Mars 压力单位必须是 Pa。**
   MACDA 的 `psurf`、`lev` 要确认：如果垂直坐标是 sigma，常见处理是 `PRES = psurf * lev`。不能用 hPa。
3. **Mars 没有地球意义上的 PMSL。**
   `real.exe` 地球路径习惯用 `PMSL` 反推 `PSFC`，但 Mars 应优先直接提供 `PSFC`，避免海平面气压逻辑污染结果。
4. **`ALBEDO12M` 在地球 WPS 里按百分比处理。**
   `real.exe` 会把 `ALBEDO12M / 100` 变成 `ALBBCK`。如果 Mars 反照率已经是 `0-1`，直接塞进去会被除以 100，结果错两个数量级。
5. **当前 `ideal.exe` 直接读 `Data/`，不经过 geogrid。**
   例如 MOLA 地形、TES albedo、thermal inertia 都是在 `module_setup_surface.F` / `module_planetary_terrain.F` 里直接读本地二进制文件。真实资料路径要么复刻到 geogrid，要么修改 `real.exe` 初始化。
6. **Mars 时间不是普通公历习惯。**
   namelist 用 `start_year/start_day`，输出时间类似 `0003-00001_00:00:00`。MACDA 的 `MY/SOY/Ls` 到 WRF 时间轴要单独处理。
7. **C-grid stagger 不能混。**
   `UU` 对应 U stagger，`VV` 对应 V stagger；`metgrid`/`real` 会处理插值，但 WPS intermediate 的变量名和 `METGRID.TBL` 标记要匹配。


# MACDA 修改清单
| 类别       | 文件/目录                                                                                               | 动作           | 用途                                                                |
| ------------ | --------------------------------------------------------------------------------------------------------- | ---------------- | --------------------------------------------------------------------- |
| MACDA 目录 | `sample/MACDA-v2/`                                                                                  | 保留并扩展     | MACDA 专用 profile 根目录                                           |
| 配置       | `sample/MACDA-v2/config.MACDA-v2.ini`                                                               | 新建           | MACDA 输入路径、输出路径、时间范围、Mars 半径、目标 pressure levels |
| 元数据表   | `sample/MACDA-v2/db/macda_meta.csv`                                                                 | 新建           | 替代`cmip6_meta.csv`，声明 MACDA 变量组和主时间频率             |
| 变量表     | `sample/MACDA-v2/db/MACDA_ATM.csv`                                                                  | 新建           | `temp/uwind/vwind/geop`等 3D 大气变量映射                       |
| 变量表     | `sample/MACDA-v2/db/MACDA_SFC.csv`                                                                  | 新建           | `psurf/tsurf/co2ice/coldust`等 2D 地表/柱变量映射               |
| 变量表     | [db/MACDA\_LEV.csv (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/db/MACDA\_LEV.csv:1)                       | 后续删除或废弃 | 当前内容是占位且有重复表头，不能作为正式表                          |
| 说明文档   | `sample/MACDA-v2/README.md`                                                                         | 新建           | 记录 MACDA v2.0 变量、单位、时间、垂直坐标和处理假设                |
| 样例信息   | `sample/MACDA-v2/macda-mro-mcs-reanalysis-head`                                                     | 保留           | 作为 MACDA 文件头信息参考                                           |
| WPS 样例   | `sample/MACDA-v2/namelist.wps`                                                                      | 修改           | 先作为样例，不直接保证可跑；需修正中文标点和`fg_name`           |
| 代码入口   | [run\_c2w.py (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/run\_c2w.py:1)                                   | 修改           | 增加`--config`参数，允许直接指定 MACDA profile 配置             |
| 适配器     | `lib/adapters/macda.py`                                                                             | 新建           | 负责打开 MACDA 单文件/多文件、Mars 时间、lat/lon 排序               |
| 适配器注册 | [lib/adapters/\_\_init\_\_.py (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/lib/adapters/\_\_init\_\_.py:1) | 修改           | `model_name=MACDA-v2`时选择`MacdaAdapter`                   |
| 主处理类   | [lib/cmip\_handler.py (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/lib/cmip\_handler.py:1)                 | 修改           | 支持 profile 内部 db、MACDA 派生变量、Mars 时间序列                 |
| 网格工具   | [utils/grid.py (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/utils/grid.py:1)                               | 小改           | 支持 Mars 半径、MACDA 原生 5 度网格                                 |
| 写文件工具 | [utils/utils.py (line 1)](/home/zy/WRF/cmip6-to-wrfinterm/utils/utils.py:1)                             | 修改           | `MAP_SOURCE`不再硬编码为`CMIP6`，改成配置项                 |
| 垂直插值   | `utils/mars_vertical.py`                                                                            | 新建           | `p = psurf * lev`，sigma 层插值到固定 Mars pressure levels      |
| 时间工具   | `utils/mars_time.py`                                                                                | 新建           | 解析`Mars_date`，生成 WRF intermediate 的`HDATE`字符串      |
| 测试       | `sample/MACDA-v2/test_macda_smoke.py`或`tests/test_macda.py`                                    | 新建           | 用 1-2 个时次生成 intermediate，检查字段和维度                      |

**MACDA 变量表建议**
| MACDA 变量 | WPS/WRF intermediate 字段 | 处理 |
|---|---|---|
| `temp(time,lev,lat,lon)` | `TT` | K，sigma -> fixed pressure |
| `uwind(time,lev,lat,lon)` | `UU` | m s-1，sigma -> fixed pressure |
| `vwind(time,lev,lat,lon)` | `VV` | m s-1，sigma -> fixed pressure |
| `geop(time,lev,lat,lon)` | `GHT` | `geop / 3.727`，m |
| `psurf(time,lat,lon)` | `PSFC` | Pa，直接写 |
| `tsurf(time,lat,lon)` | `SKINTEMP` | K，直接写 |
| 常数 0 | `QV` 或 `SPECHUMD` | 第一版可设 0 |
| `coldust` | `TAU_OD2D` | 可选，建议 `coldust / psurf * 700 Pa` |
| `co2ice` | `CO2ICE` | 可选，后续验证 `real.exe` 是否接收 |
| `dustmmr` | 暂缓 | 不能直接等同 `TRC01/TRC02`，需要物理方案确认 |
| `swflux/lwflux` | 暂缓 | 这是诊断通量，不是最小初始场 |
