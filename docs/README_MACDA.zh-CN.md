# MACDA v2.0 到 MarsWRF intermediate

MACDA 转换路径与保留下来的地球/CMIP 转换器相互独立。在本目录安装少量依赖并运行：

```bash
pip install -r requirements.macda.txt
python run_macda2w.py --dry-run
python run_macda2w.py --time-index 0
python run_macda2w.py --max-times 2
```

默认输入是 `/home/zy/WRF/DATA/macda-mro-mcs` 下的 MY28/SOY507-SOY537 文件，结果写入 `output/`。

## 时间转换

MACDA 提供三套相互独立的时间描述：以 sol 为单位的连续 `time`、`MY_Ls/Ls`，以及 `Mars_date` sol 日历。已检查的 MarsWRF 文件使用 `YYYY-DDDDD_HH:MM:SS`，其中 `DDDDD` 是从 1 开始的年内 sol 序号。转换器按照文档规定的 MACDA 五年周期 `669,668,669,668,669` 和月份长度计算该序号，再用连续 `time` 进行交叉校验。

```text
+0028-10-07T02:00:00A -> 0028-00507_02:00:00.0000
```

旧的直接字符替换只能得到 `0028-10-07_...`，这不是合法的 MarsWRF 行星日期。`XFCST` 表示从本次所选首条记录起算的火星小时数。

## 字段处理流程

当前唯一生效的变量表是 `db/MACDA-v2_CORE.csv`：

```text
temp    -> TT          sigma 层 -> 固定压力层
uwind   -> UU          sigma 层 -> 固定压力层，相对地理坐标的风
vwind   -> VV          sigma 层 -> 固定压力层，相对地理坐标的风
derived -> PRESSURE    由 metgrid 派生最终 PRES
geop    -> GHT         位势 / 火星重力加速度
psurf   -> PSFC
tsurf   -> SKINTEMP
zero    -> SPECHUMD, QV
coldust -> TAU_OD2D    归一化到 700 Pa，可选
co2ice  -> CO2ICE      可选
```

气压按 `p(k,j,i) = psurf(j,i) * lev(k)` 计算，垂直方向采用对数气压线性插值。每个写出的文件都会立刻回读，检查 HDATE、必需字段名称、维度、有限值、字节序和风场标志。

`omega`、`swflux` 和 `lwflux` 是诊断量，不是 `real.exe` 初始化状态。`dustmmr` 不能在没有科学依据的粒径分配方案时直接分给 MarsWRF 的两个沙尘粒径档，因此没有把这些变量错误标成 WRF 状态量。

已经合并完成的运行表是 `sample/metgrid/METGRID.TBL`。该表保留了 `PRESSURE -> PRES`、`TT/UU/VV/GHT/PSFC/SKINTEMP/SPECHUMD/QV` 规则，并加入 `TAU_OD2D` 和 `CO2ICE`。运行时需要把它放到 `opt_metgrid_tbl_path` 指定的目录中，文件名必须精确为 `METGRID.TBL`。
