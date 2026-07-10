# 0001 初始 MACDA 转换器脚手架

## 新增

- `.gitignore`
- `output/.gitkeep`
- `run_macda2w.py`
- `macda2wrf/__init__.py`
- `macda2wrf/config.py`
- `macda2wrf/grid.py`
- `macda2wrf/macda_reader.py`
- `macda2wrf/vertical.py`
- `macda2wrf/wrf_intermediate.py`
- `macda2wrf/converter.py`
- `conf/config.MACDA-v2.ini`
- `db/MACDA-v2_CORE.csv`
- `requirements.macda.txt`
- `README_MACDA.md`
- `change-records/README.md`
- `change-records/0001-initial-macda-scaffold.md`

## 修改

这一步没有修改已有源文件。

## 删除

这一步没有删除文件。

## 设计原因

保留复制来的 CMIP6 转换项目，不改变原有地球数据转换路径；另建入口和 `macda2wrf` 包隔离 MACDA 实现。这样设计是因为 MACDA 与 CMIP6 的数据组织和物理含义不同：MACDA 使用单个 NetCDF4/HDF5 文件、sigma 层、火星日期和火星气压量级，而 CMIP6 通常按变量拆分文件，使用地球日期和地球气压层。

第一版变量表覆盖 MarsWRF real-data 流程需要的最小字段：

```text
TT, UU, VV, PRES, GHT, PSFC, SKINTEMP, SPECHUMD, QV
```

同时输出两个可选火星字段：

```text
TAU_OD2D, CO2ICE
```

可选字段要真正进入 `real.exe` 或 `wrf.exe`，下游 MarsWRF/WPS 可能还需要相应配置或代码支持。

## 验证

- 新入口和包模块通过 `python3 -m py_compile`。
- 临时 writer 测试在 `/tmp` 成功创建一条 WRF intermediate 记录。
- 小型合成 sigma 到压力层插值测试得到有限数值结果。
- 配置和 CSV 解析正确找到 MACDA 文件、变量表、11 个字段定义、36 x 72 目标网格和 19 个压力层。
- `python3 run_macda2w.py --dry-run` 能运行到打开真实 MACDA 文件；脚手架阶段的默认 Python 环境缺少 `h5py` 和 `netCDF4`，因此在该处停止。

## 遗留问题

- 本地 Python 环境需要提供 `h5py` 或 `netCDF4`；脚手架建立时的默认环境二者都没有。
- 当前 `METGRID.TBL` 和 `realonly` Registry 不保证传递火星专用字段。
- 当时的 `hdate_strategy=mars_date` 只是把 MACDA 火星日历字符串写入 WRF intermediate，尚需对照 MarsWRF WPS 日期解析器验证后，才能进行长时间批量转换。

