# 0004 跨火星年的 HDATE 标签

## 新增

- `macda2wrf/mars_time.py`：`WrfPlanetDate`、`MARSWRF_FIXED_YEAR_SOLS`、`HDATE_STRATEGIES`、`whole_sols_since_macda_epoch`、`relabel_fixed_year`。
- `hdate_strategy = marswrf_fixed669`：以本次运行首条记录为锚点，把记录重新标注到固定 669-sol 日历上。
- `tests/test_mars_time.py` 的 `FixedYearLabelTest`：锚点、年界、年界后偏移、不跨年等价性、非法策略名。
- `docs/跨年模拟` 的「原因」「解决」两节。

## 修改

- `macda2wrf/converter.py`：`_write_one_time` 接收本次运行首条记录的 `Mars_date` 作为锚点；标签与真实 MACDA 日期不同时打印一行提示；`--dry-run` 增加 `hdate_strategy` 和首末时次的标签。
- `macda2wrf/config.py`：加载时校验 `hdate_strategy` 属于 `HDATE_STRATEGIES`，不再等到写文件时才报错。
- `macda2wrf/mars_time.py`：`MarsDate` 的 `wrf_date`/`hdate`/`filename_stamp` 改为委托给 `WrfPlanetDate`；`sols_since_macda_epoch` 改为整数部分加日内分数。
- `config/config.MACDA-v2-MY34-MY35.ini`：改用 `marswrf_fixed669`，并放开 `max_times` 以做全量转换。
- `docs/MARSWRF_AUDIT.md`、`docs/MARSWRF_AUDIT.zh-CN.md`：在跨年限制一段后补上这个不改 Fortran 的做法。

## 删除

没有删除任何文件。

## 设计原因

metgrid 的日期递推用固定 669 sol 的火星年（`WPS/metgrid/src/module_date_pack.F:17`），而 MACDA 的 sol 日历里 MY34 只有 668 sol。跨年时转换器写出 `MACDA:0035-00001_00`，metgrid 却按 669 日历去找 `MACDA:0034-00669_00`，于是 `Couldn't open file` 加 `mandatory field TT was not found`。

统一改 WPS 和 WRF 的日历最准确，但 `WRF/external/esmf_time_f90/Meat.F90` 等处也写死 669 sol，只改 WPS 会让 `real.exe`/`wrf.exe` 在同一个年界上以同样方式失败，工作量和回归风险都最大。所以先在转换器一侧统一到 MarsWRF 的固定年长：sol 间距与 `geth_newdate`/`geth_idts` 一致，metgrid 需要的文件全部存在，且不动任何 Fortran。

代价是年界之后的标签比真实 MACDA 日期晚 1 sol（WPS 认为多出来的那一天变成 MY34 的“幻影” sol 669）。时间轴仍是严格均匀的，但和观测对比时必须扣掉这 1 sol，模式内部 Ls 偏早约 0.5°。因此该策略是显式配置项，默认仍是 `marswrf_sol`，不会在不跨年的运行里改变任何行为。

## 验证

```text
$ python -m unittest discover -s tests -t .
Ran 14 tests in 0.021s
OK
```

MY34SOY665_MY35SOY027 文件年界处的 4 个时次：

```text
[macda2w] marswrf_fixed669 relabel 0035-00001_00 -> 0034-00669_00
[macda2w] marswrf_fixed669 relabel 0035-00001_22 -> 0034-00669_22
[macda2w] marswrf_fixed669 relabel 0035-00002_00 -> 0035-00001_00
output/MY34-MY35/MACDA:0034-00668_22
output/MY34-MY35/MACDA:0034-00669_00
output/MY34-MY35/MACDA:0034-00669_22
output/MY34-MY35/MACDA:0035-00001_00
```

`--dry-run` 报告 360 个时次，首末标签为 `0034-00665_02:00:00` 和 `0035-00026_00:00:00`，在 669-sol 日历下正好 359 个 7200 秒间隔。

## 遗留问题

- `namelist.wps` 的 `end_date` 需要人工改成 `0035-00026_00:00:00`；转换器不写 namelist。
- 全量 360 时次的转换和 `metgrid.exe` 尚未重跑；集群和本地输出目录里旧命名的 `MACDA:*`、`met_em.*` 必须先删掉，否则会混入不属于新序列的文件。
- 年界后的 1 sol 标签偏移要带到洞察号对比的结论里。想彻底消除，仍然要按 `docs/MARSWRF_AUDIT.zh-CN.md` 第 5 节把 WPS 和 WRF 的日历一起改成 669/668 循环，或者改用对齐好 epoch 的 MARS24/MSD 编译。
