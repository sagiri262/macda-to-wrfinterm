# 0005 跨火星年运行的稳定标签锚点

## 新增

- `macda2wrf/mars_time.py`：`MacdaRecord`、`LabelledRecord`、`make_label`、`build_label_sequence`、`check_unique_labels`、`check_label_spacing`、`check_hdate_strategy`、`label_sols`、`fixed_year_sols`、`mars_date_at`、`month_and_sol`、`whole_sols_since_epoch`、`short_mars_years_left`、`ANCHORED_HDATE_STRATEGIES`、`MarsDate.macda_string`、`WrfPlanetDate.sol_fraction`。
- `macda2wrf/converter.py`：`_anchor_index`、`_select_indices`、`_validate_indices`、`_build_labels`。
- `tests/macda_boundary_fixture.py`：用 MACDA 日历本身重建 `mro-mcs-reanalysis_mars_MY34SOY665_MY35SOY027_v2-0.nc` 的时间轴，使测试不读取那个 806 MB 的外部文件，也不依赖 `output/` 和 `public/`。
- `tests/test_mars_time.py`：`MacdaCalendarInverseTest`、`StableAnchorLabelTest`、`LabelSequenceValidationTest`，以及 `FixedYearLabelTest` 新增的两条（缺锚点、默认策略保留真实 MY35 标签）。
- `tests/test_converter.py`：`ConverterStableAnchorTest`，以及 `ConverterIndexTest` 新增的四条。
- `docs/MACDA跨年固定669兼容方案.md`。

## 修改

- `macda2wrf/mars_time.py`：`marswrf_fixed669` 现在必须显式给出锚点。`make_hdate` 不再在缺锚点时把记录当作自己的锚点，而是直接报错。`relabel_fixed_year` 的文档字符串明确它只是兼容性重贴标签。
- `macda2wrf/converter.py`：标签与 `XFCST` 的锚点改为配置 `start_index` 指定的那条记录，每次运行只读一次，不再使用当前命令的第一个 `time_index`。显式 `--time-index` 会排序去重，越界报错，早于锚点报错。整段选择在写出任何文件之前统一贴标签并校验。`--dry-run` 逐条校验并打印，不再只报首末两条。选择区间在使用真实 MACDA 标签的情况下离开一个短于 669 sol 的火星年时会给出警告。
- `config/config.MACDA-v2-MY34-MY35.ini`：只改注释。说明 `start_index = 0` 就是本次运行的稳定锚点，全量与部分重跑之间不得更改；`hdate_strategy` 处指向新文档。
- `docs/MARSWRF_AUDIT.md`、`docs/MARSWRF_AUDIT.zh-CN.md`：各改一句，把「以本次运行的首条记录为锚点」更正为配置的 `start_index`，并指向新文档。

## 删除

没有删除任何文件。

## 设计原因

记录 0004 把固定 669-sol 日历锚在「本次运行的第一条记录」上。这让同一条 NC 记录的标签取决于一条命令恰好选了哪些时次：

```text
--time-index 47                  ->  0035-00001_00
--time-index 0 --time-index 47   ->  0034-00669_00
```

于是部分重跑、断点续跑和全量运行会为同一份数据写出不同的文件名，`XFCST` 也随之漂移。锚点必须是数据集与配置的属性，而不是命令行的属性，所以现在固定为 `start_index`。单独一条记录不包含固定日历从哪里起算的信息，因此锚点是必填而不是可省：悄悄以记录自身为锚点正是这个 bug 本身，不是一个合理的退化路径。

校验从「逐文件」改为「整段序列」出于同样的理由。metgrid 是用日期算术推出文件名的，所以关键不在于每个标签格式是否正确，而在于整段序列是否唯一（没有输出文件互相覆盖）、是否严格递增、间距是否与 MACDA 的 `time` 坐标逐对一致。`build_label_sequence` 在写出第一个文件之前检查这三项，`--dry-run` 则报告整段选择。

## 验证

测试**已编写但未执行**：本次改动在「仅本地」的约束下完成，没有做任何转换、没有跑 metgrid、也没有在本地运行测试。`tests/test_mars_time.py` 与 `tests/test_converter.py` 覆盖了全量与单时次的标签一致性、MY34 sol 668 到固定 sol 669 的映射、跨年界标签连续、部分重跑时 `XFCST` 仍相对配置锚点、缺锚点报错、以及默认 `marswrf_sol` 标签不变。

本地实际执行过的静态检查：对五个改动的 Python 文件做 `ast.parse`，无超过 120 字符的行，无制表符，无「括号内单值独占一行、右括号再独占一行」的写法。

HPC 验证清单见 `docs/MACDA跨年固定669兼容方案.md` 第 7 节。

## 遗留问题

- 集群上尚未验证任何一项。在新文档第 7 节全部完成之前，本次改动只能表述为「已实现、待验证」。
- `namelist.wps` 仍写着 `end_date = '0035-00027_00:00:00'`，那是最后一条记录的真实 MACDA 日期，不是它的固定 669 标签。必须改成 `0035-00026_00:00:00`：718 h / 2 h = 359 个间隔，即 360 个时次。
- `output/MY34-MY35` 是旧的不稳定锚点下写出的输出，混了两套标签口径，不能作为验证依据。
- 年界之后的 1 sol 标签偏移仍然存在。想彻底消除，仍需按 `docs/MARSWRF_AUDIT.zh-CN.md` 第 5 节把 WPS 与 WRF 的日历一起改成 669/668 循环。
- 标识符仍沿用 `snake_case` 以与包内其余代码一致，与共享代码规范中的 lowerCamelCase 要求不符。改名属于全仓重构，超出本次改动范围。
