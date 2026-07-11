# Mars 静态地理数据对齐说明

## 1. 目标与运行文件

本次把 `/home/zy/WRF/DATA/Data` 中由 MarsWRF `ideal.exe` 直接读取的静态表面数据转换为 WPS `geogrid` 可读取的二进制 tile。原始数据没有改动。

运行时使用：

- `sample/namelist.wps.copy`：修改后的 WPS 配置，路径均为相对路径。
- `sample/geogrid/GEOGRID.TBL`：`geogrid.exe` 实际查找的精确文件名。
- `sample/geogrid/GEOGRID.TBL.ARW.copy`：保留 `.copy` 后缀的修改结果。
- `sample/WPS_GEOG_MARS/*/index`：实际读取的索引文件。
- `sample/WPS_GEOG_MARS/*/00001-00360.00001-00180`：实际读取的全球 tile。
- `sample/WPS_GEOG_MARS/build_manifest.json`：每个 tile 的大小、层数、值域和 SHA-256。

`index.copy` 是此前的模板，仍然保留，但运行时不会读取它。

## 2. 与当前 MarsWRF 配置的对应关系

转换脚本读取 `em_global_mars/namelist.input` 所采用的组合：

| MarsWRF 设置 | 使用的数据 | 处理方式 |
|---|---|---|
| `topo64=.true.` | `Data/topo/topo_latlon/topo64` 下四个 MOLA 象限 | 1/64 度按 64×64 平均为 1 度，并拼成 0～360°E 全球网格 |
| `alb_my=1` | `AlbMY24map.bin`，极区由旧 `tes_albedo_filled.dat` 补充 | 按 MarsWRF 相同的经向翻转及 IAU 1994→2000 偏移处理，再聚合为 1 度 |
| `ti2007=0` | `tes_inertia_filled.dat` | 保持旧 TES 热惯量选择，不能替换成 DB/NB 2007 产品 |
| Mars 表面粗糙度 | `roughness.dat` | 按源码中的 `0.15/198` 换算为米 |
| 地下冰比例、深度 | `grs_subsurface_ice.dat`、`h2odepth.dat` | 保留原值和 `-9999` 缺测标记 |

TES 数据采用 IAU 1994 经度，而 MOLA/WRF 使用 IAU 2000。转换脚本复现源码中的 `LON(1994)=LON(2000)-0.271`，最后统一插值到中心为 `0.5, 1.5, ..., 359.5°E` 的网格。

## 3. geogrid 字段与数据来源

| geogrid 字段 | WPS_GEOG_MARS 目录 | 来源或设定 |
|---|---|---|
| `HGT_M` | `mola_topography` | topo64 MOLA |
| `ALBEDO`、`ALBBCK` | `mars_albedo` | 当前 MY24 TES 反照率，单位为 0～1 |
| `ALBEDO12M` | `mars_albedo12m` | 同一 MY24 场重复 12 层，geogrid 中为百分数；`real.exe` 会除以 100 |
| `THC`、`THCBCK` | `mars_thermal_inertia` | `ti2007=0` 对应的旧 TES 场，单位 tiu |
| `ZNT` | `mars_roughness` | MOLA pulse-width 粗糙度，单位 m |
| `grd_ice_pc` | `mars_soil_ice_pc` | GRS 地下冰体积分数 |
| `grd_ice_dp` | `mars_soil_ice_dp` | 地下冰层顶深度，单位 m |
| `EMISS`、`EMBCK` | `mars_emissivity` | 与 ideal 表面初始化一致，固定为 1.0 |
| `LANDUSEF` | `mars_surface_class` | 全火星陆面为类别 1；保留 24 类轴以匹配 WRF 默认 `num_land_cat=24` |
| `SOILCTOP`、`SOILCBOT` | `mars_soil_class` | 统一类别 1；保留 16 类轴以匹配默认 `num_soil_cat=16` |
| `GREENFRAC` | `mars_greenfrac` | 无植被，12 层均为 0 |
| `SOILTEMP` | `mars_deep_soil_temp` | `Data` 中无相应气候场，兼容值为 200 K |
| `SNOALB` | `mars_snow_albedo` | `Data` 中无独立雪反照率气候场，兼容值为 0% |
| `H2OICE` | `mars_h2oice` | 当前 `mp_physics=48` 的 ideal 初始化为 0 |

表中明确写为全局固定值或兼容值的字段不是观测产品，不能解释为从 TES/MOLA 反演得到的数据。

## 4. 重新生成

在项目根目录执行：

```bash
python utils/build_mars_geogrid_data.py \
  --data-root /home/zy/WRF/DATA/Data \
  --output-root sample/WPS_GEOG_MARS
```

在 HPC 上需要把 `DATA/Data`、项目目录及生成后的 `WPS_GEOG_MARS` 一起同步，或者在 HPC 上用其实际路径重新执行脚本。脚本要求 NumPy，并会先检查各源文件的大小、大端字节序 Fortran 记录头和物理值域。

## 5. HPC 运行位置

相对路径配置要求从 `sample` 目录运行：

```bash
cd /public/home/proj_kcchow/zhaoy/marswrf/macda-to-wrfinterm/sample
ln -sf ../../WPS/geogrid.exe ./geogrid.exe
cp namelist.wps.copy namelist.wps
test -r geogrid/GEOGRID.TBL
test -r WPS_GEOG_MARS/mola_topography/index
test -r WPS_GEOG_MARS/mola_topography/00001-00360.00001-00180
./geogrid.exe
```

如果从其他运行目录启动，必须同步修改 `geog_data_path` 和 `opt_geogrid_tbl_path`。仅存在 `index.copy` 或 `GEOGRID.TBL.ARW.copy` 不够，因为 WPS 读取的是无 `.copy` 后缀的精确文件名。

本次没有编译或执行 `geogrid.exe`。二进制 tile 已按 WPS 的 x-fastest、big-endian 整数布局生成，并通过尺寸、值域、哈希及重复生成一致性检查。
