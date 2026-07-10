# 本文内容均由 ChatGPT 生成，请务必甄别内容

下面按你的 3 个问题解释。这个仓库本质是一个 **ETL 工具**：把 CMIP6 或类似 NetCDF 数据读取出来，按 WRF/WPS 需要的变量名、层次、单位、网格和二进制格式，写成 WRF intermediate file，后续给 `metgrid.exe` 使用。

**1. conf 和 db**
`conf/*.ini` 是“运行参数配置”。它回答的是：这次跑哪个模型、数据在哪里、情景是什么、输出哪个时间段、结果写到哪里。

例如 [conf/config.CESM2.ini](/home/zy/WRF/cmip6-to-wrfinterm/conf/config.CESM2.ini:1) 里有：

```ini
[INPUT]
input_root=./sample/CESM2/
model_name=CESM2
scenario=hist
esm_flag=r11i1p1f1
grid_flag=gn

[OUTPUT]
etl_strt_ts = 199001010000
etl_end_ts  = 199001020000
output_root = ./output/
```

为什么要 `.ini`：因为这些是“会经常换，但不应该写死在代码里”的东西。比如换模型、换情景、换输入目录、换输出时间段时，只改配置，不改函数。`.ini` 不是唯一选择，YAML/TOML 也可以，但 Python 标准库自带 `configparser`，所以实现简单。

`conf/logging_config.ini` 是日志配置，控制日志打印到屏幕和 `c2w.log`。

`db/*.csv` 更接近 WPS 里的 **Vtable/变量表**。它不是单纯做变量名一一映射，而是同时承担 4 件事：

1. 源变量名是什么：`src_v`，比如 `ta`
2. 写给 WRF 的变量名是什么：`aim_v`，比如 `TT`
3. 这个变量是什么类型：`3d`、`2d`、`2d-soil-mon`、`2d-fixed` 等
4. 变量在哪种垂直坐标上：`Lev` 表示模式层/混合层，`PlevPt` 表示气压层

例如 [db/CESM2_6hrLev.csv](/home/zy/WRF/cmip6-to-wrfinterm/db/CESM2_6hrLev.csv:1)：

```csv
src_v,aim_v,units,type,lvlmark,desc
ta,TT,K,3d,Lev,3-d air temperature
ua,UU,m s-1,3d,Lev,3-d wind u-component
ps,PSFC,Pa,2d,Lev,Surface pressure
```

这确实包含“变量名映射”，但不止于此。`type` 会直接决定代码怎么处理：是三维插值、二维插值、土壤层重映射、固定场，还是日/月频率变量。

还有一个核心表是 [db/cmip6_meta.csv](/home/zy/WRF/cmip6-to-wrfinterm/db/cmip6_meta.csv:1)。它决定某个模型、某个 scenario 要加载哪些变量组：

```csv
model_name,scenario,variable_group,naming_convention,var_frq,table_id
CESM2,hist,6hrLev,,6h*,6hrLev
CESM2,hist,3hr,,3h,3hr
CESM2,hist,Lmon,,1M,Lmon
```

其中 `6h*` 的星号表示主输出时间频率。也就是说，程序按这个频率生成 `out_time_series`。

**2. Python 文件职责**
工程化地看，这个仓库分成 5 层：

```text
run_c2w.py
  -> lib/cfgparser.py
  -> lib/cmip_handler.py
       -> lib/adapters/*.py
       -> utils/*.py
  -> output/WRF intermediate files
```

### [run_c2w.py](/home/zy/WRF/cmip6-to-wrfinterm/run_c2w.py:23) 是入口脚本。

`main_run()` 做这些事：

1. 读取日志配置
2. 解析命令行参数 `-m MODEL`
3. 找到 `conf/config.<MODEL>.ini`
4. 调用 `lib.cfgparser.read_cfg()` 读取配置
5. 创建 `CMIPHandler`
6. 遍历输出时间序列
7. 每个时刻调用 `parse_data()`
8. 再调用 `write_wrfinterm()` 写 WRF intermediate 文件

### [loop_mon_span.py](/home/zy/WRF/cmip6-to-wrfinterm/loop_mon_span.py:1) 
BCMM 的批处理脚本。它按 1 到 12 月循环修改 `conf/config.BCMM.ini` 的起止时间，然后调用 `python3 run_c2w.py -m BCMM`。这不是核心框架，只是自动跑多个月的辅助脚本。

### [lib/cfgparser.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/cfgparser.py:6)

- `read_cfg(config_file)`：用 `configparser` 读 `.ini`
- `write_cfg(cfg_hdl, config_fn)`：把配置写回文件

### [lib/cmip_handler.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/cmip_handler.py:30) 
该程序是核心程序。`CMIPHandler` 负责“从数据到 WRF intermediate”的主流程。

`CMIPHandler.__init__()` 初始化整套转换任务：读输入/输出配置、确定模型名、scenario、ETL 时间段、输出目录、输出网格、adapter、meta 表、加载数据、准备 WRF intermediate 模板。

`_load_vtable(group_name)` 读取并缓存 `db/<MODEL>_<GROUP>.csv`。

`_build_meta(in_cfg)` 从 `db/cmip6_meta.csv` 里筛选当前模型和 scenario 的变量组，并根据带 `*` 的频率生成输出时间序列。

`_extract_hybrid_coeffs(ds)` 从混合 sigma 层数据里提取 `ap`、`b`、`ps`，用于把模式层转换到标准气压层。

`_load_sftlf_helper()` 尝试加载 `sftlf_*.nc` 陆地比例掩膜，用于填补陆地/海洋变量的缺测值。

`_emit_namelist_hints()` 如果某模型选择跳过土壤变量，则输出一个 JSON 提示用户如何改 WRF `namelist.input`。

`_load_cmip_data()` 根据 meta 表和 Vtable 打开所有需要的 NetCDF 变量，并做基础整理：时间裁剪、变量去重、`lev`/`plev` 坐标统一、CESM2 层顺序修正、土壤深度边界读取等。

`_interp_to_grid(da, with_plev=False)` 把变量插值到统一输出网格，必要时也插值到统一气压层。

`_fill_2d(da, src_kind='full')` 用最近邻填补二维变量中的 NaN，尤其处理海洋变量、陆地变量跨海岸线的问题。

`parse_data(tf)` 是每个时刻的数据转换核心：选取时间片、做单位转换、三维变量转气压层、二维变量填缺测、土壤变量重映射、最后把结果放进 `self.outfrm`。

`write_wrfinterm(tf, tgt)` 把 `self.outfrm` 写成 WRF intermediate 二进制文件。`tgt='main'` 写大气/陆面主文件，`tgt='sst'` 写 SST 文件。

[lib/adapters/_base.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/adapters/_base.py:19) 定义抽象基类 `ModelAdapter`。它把“文件在哪里、怎么打开、时间坐标怎么选”从主处理逻辑里拆出去。

`ModelAdapter` 的方法：

- `__init__()` 保存模型名、输入根目录、scenario、ETL 时间
- `_files_for()` 抽象方法，子类必须实现，用来找文件
- `open_for()` 根据变量打开对应 Dataset，并缓存
- `time_to_index()` 把普通 Python 时间转换成 xarray 可用的时间选择值；no-leap 日历会转成 `cftime.DatetimeNoLeap`
- `close()` 关闭缓存的 Dataset
- `_cache_key()` 决定缓存粒度
- `_open()` 用 `xarray.open_dataset/open_mfdataset` 打开文件

[lib/adapters/cmip6.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/adapters/cmip6.py:28) 定义 `Cmip6Adapter`，负责标准 CMIP6 文件命名：

```text
<var>_<table>_<model>_<scenario>_<member>_<grid>_<time-range>.nc
```

它支持两种发现方式：

- `exact`：配置里有 `cmip_strt_ts/cmip_end_ts`，按精确文件名找
- `glob`：没有起止后缀时，用通配符找，更适合真实 CMIP6 数据文件分块不统一的情况

`_file_scenario()` 把配置里的 `hist` 转成文件名里的 `historical`。

`_table_name()` 从 `cmip6_meta.csv` 的 `table_id` 找 CMIP6 表名，比如 `6hrLev`、`Amon`、`Lmon`。

[lib/adapters/bcmm.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/adapters/bcmm.py:15) 定义 `BcmmAdapter`。BCMM 不是标准“一变量一文件”的 CMIP6 布局，而是一个月一个文件、一个文件里有很多变量。所以它设置：

```python
one_ds_per_group = True
soil_packed_4d = True
```

`_files_for()` 根据 `cmip6_meta.csv` 里的 `naming_convention` 生成文件名。

[lib/adapters/__init__.py](/home/zy/WRF/cmip6-to-wrfinterm/lib/adapters/__init__.py:43) 主要提供工厂函数 `make_adapter()`。它根据模型名选择 `BcmmAdapter` 或 `Cmip6Adapter`，并处理 `calendar`/`cftime`。

[utils/grid.py](/home/zy/WRF/cmip6-to-wrfinterm/utils/grid.py:18) 定义 `OutputGrid`，表示输出网格。默认是全球 1 度网格、14 个标准气压层、4 个 WRF 土壤层。`from_config()` 允许通过 `.ini` 的 `[OUTPUT]` 改 `lat/lon/nlat/nlon/plev_hpa/soil_layers`。

[utils/soil.py](/home/zy/WRF/cmip6-to-wrfinterm/utils/soil.py:44) 负责土壤层处理：

- `parse_wrf_soil_label(label)`：把 `ST040100` 解析成 40-100 cm
- `remap_soil_layer()`：把源数据的土壤深度层按重叠厚度加权平均到 WRF 四层
- `read_depth_bnds()`：从 NetCDF 的 CF convention 里读取 `depth_bnds`

[utils/utils.py](/home/zy/WRF/cmip6-to-wrfinterm/utils/utils.py:26) 是通用工具：

- `throw_error()`：记录错误并退出
- `write_log()`：写日志
- `gen_wrf_mid_template()`：生成 WRF intermediate 单条记录的模板
- `write_record()`：按 WRF intermediate 格式写 Fortran unformatted record
- `fill_nan_2d_nearest()`：二维最近邻填 NaN
- `hybrid2pressure()`：把混合 sigma 层变量转换到标准气压层

这几个文件组合起来实现的功能是：

```text
配置 .ini
  -> 选择模型、时间、目录、输出网格
db/cmip6_meta.csv
  -> 选择变量组、频率、CMIP6 table
db/<MODEL>_<GROUP>.csv
  -> 选择变量、WRF 字段名、变量类型、垂直层语义
adapter
  -> 找文件、开文件、处理日历
CMIPHandler
  -> 读数据、插值、填缺测、转层、处理土壤、写文件
utils
  -> 网格、土壤、二进制 WRF intermediate 写入
```

**3. 迁移性判断**
这个仓库有一定迁移性，但不是“任何新数据只改几个值就能跑”。

如果你的新再分析数据满足这些条件，通常是小改：

- 文件可以整理成类似 CMIP6 的命名方式，或者只需要新建一个简单 adapter
- 变量坐标叫 `time/lat/lon/plev` 或 `time/lev/lat/lon`
- 三维变量已经在气压层，或有 `ap/b/ps` 可转气压层
- 变量只需要简单单位转换
- WRF 需要的变量都能在 Vtable 里对应上
- 只想改输出气压层，可以在 `[OUTPUT]` 加 `plev_hpa=1000,925,850,...`
- 只想改输出网格，可以在 `[OUTPUT]` 加 `lat_start/lat_end/nlat/lon_start/lon_end/nlon`

这类情况主要改：

1. 新建 `conf/config.<NEW>.ini`
2. 在 `db/cmip6_meta.csv` 加新模型行
3. 新建 `db/<NEW>_<GROUP>.csv`
4. 必要时加 `calendar=noleap` 或类似配置

如果你的新数据不符合这些假设，就可能要改代码：

- 文件命名完全不同：需要新写 adapter，类似 `BcmmAdapter`
- 坐标名不是 `lat/lon/time/plev/lev`：要在 adapter 或 handler 里标准化
- 是 curvilinear grid、rotated pole、非规则网格：当前 `xarray.interp(lat, lon)` 不够，需要重写水平插值
- 垂直坐标不是气压层也不是 CMIP6 hybrid `ap/b/ps`：需要新增垂直转换逻辑
- 单位转换复杂，例如湿度是相对湿度但 WRF 需要比湿：需要加物理转换函数
- 土壤层定义和 WRF 四层差异很大，且没有 `depth_bnds`：需要重新判断土壤处理策略
- 变量缺失，需要诊断计算，比如没有 `zg/GHT`、没有 `psl/PMSL`、没有近地面风：需要新增计算方法或接受 WRF fallback

所以结论是：**对标准 CMIP6 模型或“伪装成 CMIP6 结构”的数据，迁移成本较低；对任意新再分析资料，迁移性取决于文件布局、坐标系统和变量物理含义。**

你提到“更换 plev 等参数”：`plev_hpa` 这类输出层可以算小改；变量名映射也算小改；但如果源数据的垂直坐标、网格结构、文件组织方式不同，就不是只改参数，需要动 adapter 或 `CMIPHandler.parse_data()` 的处理逻辑。