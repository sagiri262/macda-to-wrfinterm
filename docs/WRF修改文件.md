## WRFV3 修改文件

| Fortran 文件 | 使用程序 | 修改内容 | 修改后编译 |
|---|---|---|---|
| [main/real_em.F](/home/zy/WRF/marswrf/WRFV3/main/real_em.F:567) | `real.exe` | 将日期构造分支改为 `WRF_PLANET`，生成 `YYYY-DDDDD_HH:MM:SS` 火星日期 | `./compile em_real` |
| [share/input_wrf.F](/home/zy/WRF/marswrf/WRFV3/share/input_wrf.F:178) | `real.exe`、`wrf.exe` 共用 | 按火星格式读取 `SIMULATION_START_DATE`；用 `start_year + start_day` 构造 `first_date_nml`，使其与 `met_em` 日期匹配 | `./compile em_real` |
| [dyn_em/module_initialize_real.F](/home/zy/WRF/marswrf/WRFV3/dyn_em/module_initialize_real.F:880) | `real.exe` | 避免使用地球公历月份处理火星 `GREENFRAC/ALBEDO12M`；增加 `WRF_PLANET` 年度气候场处理分支 | `./compile em_real` |

这三个是本次接入真正修改的 WRF Fortran 源文件。应修改大写扩展名的 .F 源文件，不要修改编译生成的 .f90。`./compile em_real` 会重新生成/链接 real.exe 和对应的 wrf.exe，不需要执行 clean。`share/module_date_time.F`、`module_mars24.F`、`module_planet_utilities.F`、`module_model_constants.F` 和 `module_setup_v_grid.F` 是 MarsWRF 原有基础设施，本次主要是参考和调用，没有列入修改文件。

## WPS 修改文件
| Fortran 文件 | 使用程序 | 修改内容 | 修改后编译 |
|---|---|---|---|
| [geogrid/src/constants_module.F](/home/zy/WRF/marswrf/WPS/geogrid/src/constants_module.F:26) | `geogrid.exe` | 定义 MarsWRF 火星半径 `MARS_RADIUS_M=3389.92e03` 和周长 | `./compile geogrid` |
| [geogrid/src/module_map_utils.F](/home/zy/WRF/marswrf/WPS/geogrid/src/module_map_utils.F:231) | `geogrid.exe` | 投影结构的 `re_m` 改用火星半径 | `./compile geogrid` |
| [geogrid/src/gridinfo_module.F](/home/zy/WRF/marswrf/WPS/geogrid/src/gridinfo_module.F:427) | `geogrid.exe` | 全球及区域经纬网格的 `dx/dy` 距离计算改用火星半径 | `./compile geogrid` |
| [metgrid/src/constants_module.F](/home/zy/WRF/marswrf/WPS/metgrid/src/constants_module.F:26) | `metgrid.exe` | 为 metgrid 单独定义相同的火星半径和周长 | `./compile metgrid` |
| [metgrid/src/module_map_utils.F](/home/zy/WRF/marswrf/WPS/metgrid/src/module_map_utils.F:231) | `metgrid.exe` | metgrid 投影计算中的 `re_m` 改用火星半径 | `./compile metgrid` |
| [metgrid/src/read_met_module.F](/home/zy/WRF/marswrf/WPS/metgrid/src/read_met_module.F:143) | `metgrid.exe` | 读取中间数据时，将网格元数据中的 `earth_radius` 设置为火星半径（单位 km） | `./compile metgrid` |
| [metgrid/src/module_date_pack.F](/home/zy/WRF/marswrf/WPS/metgrid/src/module_date_pack.F:17) | `metgrid.exe` | 日期解析和加减改为 `YYYY-DDDDD_HH:MM:SS`；使用固定 `669 sol` 火星年 | `./compile metgrid` |