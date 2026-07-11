# WPS geogrid 文件盘点

盘点日期：2026-07-11

盘点范围：`/home/zy/WRF/marswrf/macda-to-wrfinterm` 以及 `/home/zy/WRF` 中可复用的 WPS、MarsWRF 和静态数据文件。

## 1. 目标目录原始状态

| 项目 | 原始状态 | 结论 |
|---|---|---|
| `sample/namelist.wps` | 存在 | 含全角标点、缺逗号、错误前缀、地球 Lambert 参数和不完整日期字段，不能作为有效 namelist |
| `sample/geogrid/GEOGRID.TBL.ARW` | 存在 | 与 `marswrf/WPS/geogrid/GEOGRID.TBL.ARW` 完全相同，主要是旧版地球静态字段表 |
| 根层 `geogrid/GEOGRID.TBL.ARW` | 不存在 | 新增 `.copy` 结果，但不创建或覆盖无后缀原文件 |
| Mars 静态数据 `index` | 不存在 | 新增 11 份 `index.copy` 模板；没有伪造 tile 数据 |
| `geogrid.exe` | 不存在 | 从 WPS 4.6.0 复制参考二进制为 `geogrid.exe.copy` |
| NetCDF I/O 动态库 | 项目内不存在，系统中存在 | 把二进制直接依赖的 `libnetcdff` 和 `libnetcdf` 复制为留痕副本 |

## 2. GEOGRID.TBL 来源核对

样例原表并不是 `/home/zy/WRF/WPS-4.6.0/geogrid/GEOGRID.TBL.ARW` 的副本：

```text
sample 原表 SHA256       377a97061e1eab284498c5de2b3c8793d495dc75cbb0166246de8914c905a39e
marswrf/WPS 原表 SHA256  377a97061e1eab284498c5de2b3c8793d495dc75cbb0166246de8914c905a39e
WPS 4.6.0 官方表 SHA256  dd1237f317ed5162d03cd289516e8bbd9e88375f053eb7093e26456a20b10a33
```

因此，样例原表实际来自旧 `marswrf/WPS`。真正的 WPS 4.6.0 表已另存为：

```text
sample/geogrid/GEOGRID.TBL.WPS-4.6.0.reference.copy
```

## 3. 找到的 Mars 原始静态资料

工作区已有可作为转换源的数据，但它们不是 WPS geogrid binary tile，不能直接通过增加 `index` 使用：

| 目标字段 | 已有原始资料 | 当前格式 |
|---|---|---|
| `HGT_M` | `marswrf/em_global_mars/Data/topo/topo_latlon/topo64/*.img` | 四块 MOLA 1/64 度原始影像 |
| `HGT_M` 备用 | `marswrf/em_global_mars/Data/old_surface/megt90n000fb.img` | MOLA 全球原始影像 |
| `ALBEDO/ALBBCK` | `marswrf/em_global_mars/Data/albedo/AlbMY24map.bin` 等 | MarsWRF 自有二进制布局 |
| `THC/THCBCK` | `marswrf/em_global_mars/Data/thermal_inertia/DBmap2007.bin`、`NBmap2007.bin` | MarsWRF 自有二进制布局 |
| `EMISS/EMBCK` | 未找到独立全球 geogrid 产品 | 需要观测产品或有依据的常数场 |
| `H2OICE` | `marswrf/em_global_mars/Data/subsurface/` 和 polar 数据 | 不是直接可用的 geogrid 全球 tile |

必须先确认字节序、维度、缺测值、经度起点、纬度方向和单位，再转换成 WPS tile，例如 `00001-00360.00001-00180`。

## 4. 可执行文件和 I/O 库

找到的唯一已编译 `geogrid.exe` 为：

```text
/home/zy/WRF/WPS-4.6.0/geogrid/src/geogrid.exe
SHA256 3d27cdbc81d284fe0e6211c57131bd2f5ac50c6583312b3886d486d90751eb68
```

其直接 ELF 依赖为 `libnetcdff.so.7`、`libgfortran.so.5`、`libm.so.6` 和 `libc.so.6`；`libnetcdff` 继续依赖 `libnetcdf`、HDF5、curl 等系统库。`ldd` 检查时所有依赖都能解析。

项目内留痕副本：

```text
sample/geogrid/geogrid.exe.copy
sample/geogrid/lib/libnetcdff.so.7.1.0.copy
sample/geogrid/lib/libnetcdf.so.19.copy
```

这些 `.copy` 名称不能直接满足 ELF SONAME，只用于保留来源。实际验证时可以继续使用系统中已经可解析的原库。

## 5. 两个源码级阻塞点

1. WPS 4.6.0 `geogrid.exe` 的 `constants_module` 和 `module_map_utils` 使用地球半径 `6370000 m`。全球 5 度经纬坐标仍可生成，但输出的米制 `DX/DY` 和地图尺度不是火星值。
2. 旧 `marswrf/WPS` 的日期模块按 `YYYY-MM-DD` 解析，而目标 intermediate 和 MarsWRF 使用 `YYYY-DDDDD`。本次 namelist 采用正确的火星目标日期，但 WPS 日期源码仍需另行修改和编译。

本任务按要求不编译、不运行，因此以上问题被明确记录，没有通过地球日期或地球半径伪装规避。
