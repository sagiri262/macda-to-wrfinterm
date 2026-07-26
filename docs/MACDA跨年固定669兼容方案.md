# MACDA 跨火星年的固定 669 sol 兼容方案

本文只讨论一件事：当 MACDA v2.0 强迫场跨越一个 **668 sol** 的火星年边界时，如何让
使用 **固定 669 sol** 年长的 MarsWRF / WPS 仍然能找到并读入中间文件。

适用范围：`config/config.MACDA-v2-MY34-MY35.ini`（MY 34 SOY 665 → MY 35 SOY 027）
以及此后任何跨年的 MACDA 转换。默认的单年配置 `config/config.MACDA-v2.ini` 不受
本方案影响。

---

## 1. 冲突从哪里来

### 1.1 MACDA / Montabone 的真实 sol 日历

MACDA v2.0 的 `Mars_date` 使用 Montabone 等人的 sol 日历，年长度以五年为周期
循环 `669, 668, 669, 668, 669`：

```
mars_year_length(year) = 668 if year % 5 in {2, 4} else 669
```

因此 **MY 34 只有 668 sol**，MY 35 有 669 sol。真实的年边界是：

| 真实 MACDA 日期 | MarsWRF 形式 |
| --- | --- |
| `+0034-12-56T22:00:00A` | `0034-00668_22:00:00` |
| `+0035-01-01T00:00:00A` | `0035-00001_00:00:00` |

也就是说，**真实日历里根本不存在 `0034-00669`**。

### 1.2 传统 MarsWRF / WPS 的固定年长度

WPS 与 WRF 的火星日期算术把一年固定为 669 sol、一 sol 固定为 86400 s：

- `WPS/metgrid/src/module_date_pack.F:17` → `INTEGER, PARAMETER :: MARS_YEAR_SOLS = 669`
- `WPS/util/src/module_date_pack.F` 同一常量
- `WRF/external/esmf_time_f90/Meat.F90` 同一约定

metgrid **不会**去 glob 目录里已有的文件名。它用 `namelist.wps` 的
`start_date`、`interval_seconds`，以 `geth_newdate(start_date, t * interval_seconds)`
逐步推算出第 `t` 个时次应有的文件名，再用 `geth_idts` 算出总时次数。推算全部走
固定 669 sol 年。

### 1.3 结果

真实标签下，`0034-00668_22:00:00` 的下一个时次是 `0035-00001_00:00:00`；
但 metgrid 在固定 669 sol 日历上前进 2 h，得到的是 `0034-00669_00:00:00`。
于是出现：

```
ERROR: Couldn't open file MACDA:0034-00669_00 for input.
```

---

## 2. 本方案的定位：兼容性重贴标签，不是日历修正

`hdate_strategy = marswrf_fixed669` 的作用是：**把真实 MACDA 记录重新贴上固定
669 sol 日历的标签**，使 metgrid 推算出的文件名与磁盘上实际存在的文件名一一对应。

必须明确：

- 它 **不修正** MACDA 日历，也 **不修正** WPS 的固定 669 sol 假设。
- 它 **不在时间上移动数据**。每条记录的物理内容、`time` 坐标、Ls 都没有变，
  变的只是写在 HDATE 与文件名里的那串标签。
- 真实日期始终保留在 NC 文件里，转换器的 dry-run 会把
  `真实 Mars_date -> 真实 MarsWRF 日期 -> 输出标签` 三者同时打印出来。
- 它只是让 WPS 能跑起来的工程手段。真正的修正应当是把 WPS/WRF 的
  `MARS_YEAR_SOLS` 换成随年份变化的 MACDA 年长度，那是另一项独立工作。

---

## 3. 稳定锚点的定义

### 3.1 为什么必须有锚点

固定 669 sol 标签只能相对某个参考记录定义：单独看一条 MACDA 记录，无法知道固定
日历应该从哪里开始数。

历史缺陷正在这里：转换器曾用**当前命令选中的第一个 `time_index`** 当锚点，于是

- `--time-index 47` 单独重跑 → 标签 `0035-00001_00`
- `--time-index 0 --time-index 47` → 同一条记录标签 `0034-00669_00`

同一条 NC 记录写出两个不同的文件名，断点续跑与部分重跑都不可信。

### 3.2 现在的定义

> **锚点 = 配置文件 `start_index` 所指的那一条记录。**

- 锚点由 `config/config.MACDA-v2-MY34-MY35.ini` 的 `start_index = 0` 显式给出，
  在整个运行期间只读取一次，与本次命令选了哪些 `--time-index` 无关。
- 锚点记录保留它自己的真实 MACDA 标签。
- 其后每过 1 个真实 sol，标签在固定 669 sol 日历上前进 1 sol：

```
delta_sols  = whole_sols_since_macda_epoch(record) - whole_sols_since_macda_epoch(anchor)
anchor_sols = anchor.year * 669 + anchor.sol_of_year - 1
index       = anchor_sols + delta_sols
year, sol   = divmod(index, 669)          # sol_of_year = sol + 1
```

- `XFCST` 同样相对这一个锚点计算：`(time[i] - time[start_index]) * 24` 火星小时。
- `marswrf_fixed669` 缺少锚点时 **直接报错**，绝不退化为“以自己为锚点”。
- 显式给出的 `--time-index` 若早于锚点，转换器直接报错（否则 `XFCST` 为负）。
  要处理更早的记录，应当下调配置里的 `start_index`，而不是在命令行绕过它。

**因此 `start_index` 一旦确定就不可再改。** 改动它等于改动全部输出文件名。

---

## 4. 真实日期与输出标签对照

配置：`start_index = 0`，锚点 = `+0034-12-53T02:00:00A`（`0034-00665_02`），
共 360 条记录，间隔 2 火星小时。

| time_index | 真实 MACDA 日期 | 真实 MarsWRF 日期 | 输出标签（文件名戳） | XFCST (h) | 说明 |
| --- | --- | --- | --- | --- | --- |
| 0 | `+0034-12-53T02:00:00A` | `0034-00665_02` | `0034-00665_02` | 0 | 锚点，标签不变 |
| 11 | `+0034-12-54T00:00:00A` | `0034-00666_00` | `0034-00666_00` | 22 | 边界前，标签不变 |
| 46 | `+0034-12-56T22:00:00A` | `0034-00668_22` | `0034-00668_22` | 92 | MY 34 最后一条 |
| 47 | `+0035-01-01T00:00:00A` | `0035-00001_00` | `0034-00669_00` | 94 | **幻影 sol 669** |
| 58 | `+0035-01-01T22:00:00A` | `0035-00001_22` | `0034-00669_22` | 116 | 幻影 sol 的最后一条 |
| 59 | `+0035-01-02T00:00:00A` | `0035-00002_00` | `0035-00001_00` | 118 | 边界后偏移 1 sol |
| 359 | `+0035-01-27T00:00:00A` | `0035-00027_00` | `0035-00026_00` | 718 | 最后一条 |

要点：

- 索引 0–46（真实 MY 34 段）标签与真实日期 **完全一致**。
- 索引 47 起，标签比真实 MACDA 日期 **早 1 sol**（真实 MY 35 sol *N* 被贴成
  MY 34 sol 669 或 MY 35 sol *N−1*）。这个偏移是常量，不再增长。
- 标签序列在固定 669 sol 日历上严格等间距（每步 2/24 sol），
  与 NC 的 `time` 坐标步长逐对一致，这正是 metgrid 需要的。
- `XFCST` 与真实时间轴一致，永远是 `2 h × time_index`，不受重贴标签影响。

---

## 5. 边界后 1 sol 偏移的科学解释限制

这是本方案唯一的实质代价，使用结果时必须记住：

1. **模拟时间轴本身没有错。** 相邻时次间隔恒为 7200 s，数据顺序、时长、
   `XFCST` 都正确，动力学积分不受影响。
2. **边界之后的日期标签不能当作真实火星日期引用。** `met_em` / `wrfout` 上写着
   `0035-00001_00` 的那一场，物理上其实是真实 MY 35 **sol 2** 00:00 MTC。
3. **由标签反推的 Ls 会偏约 0.5°。** 1 sol ≈ 0.5° Ls（随季节变化）。任何按 Ls
   分箱、与观测配对、或与其他数据集交叉比较的分析，都必须换回真实日期，
   不能直接用文件名或 HDATE。
4. **不要跨年边界做“同一 sol_of_year 逐年对比”。** 边界两侧的标签口径不同：
   前段是真实 MACDA 日历，后段带 1 sol 偏移。
5. **WRF 内部的季节强迫按标签走。** 太阳赤纬、日照等由模式日期驱动，边界后
   相当于提前了 1 sol 的季节相位。对数天量级的短期模拟影响很小，
   但长期（数十 sol 以上）跨年模拟应当评估这一项，或改用真正的可变年长度日历。
6. **每跨一个 668 sol 的年，偏移再累加 1 sol。** 本文件只跨一次边界，偏移为 1 sol；
   更长的多年模拟会线性累积，届时固定 669 sol 方案不再适用。

若研究目标本身依赖精确的火星日期或 Ls，应放弃本兼容方案，改为修正
WPS/WRF 的年长度常量。

---

## 6. 本次改动范围（全部为本地代码改动）

本次工作 **只改本地代码与本地配置**，没有连接 HPC，没有运行 geogrid / metgrid /
real.exe / wrf.exe，没有编译任何 Fortran/C 代码，没有执行完整 NC 转换，也没有
生成新的 `output/` 或 `met_em`。

- `macda2wrf/mars_time.py` —— 锚点化标签、MACDA 日历反函数、整段序列校验。
- `macda2wrf/converter.py` —— 锚点固定为 `start_index`；索引校验；整段贴标签；
  dry-run 逐条打印与校验；`XFCST` 相对锚点。
- `config/config.MACDA-v2-MY34-MY35.ini` —— 仅补充注释，说明
  `start_index = 0` 就是锚点、`marswrf_fixed669` 只是兼容策略。
- `tests/` —— 单元测试与合成时间轴（不读取 806 MB 的外部 NC，
  不依赖 `output/` 与 `public/`）。测试 **已编写但未在本地执行**，交由 HPC 运行。

---

## 7. HPC 验证清单

以下命令需在 HPC 上执行。**在全部通过之前，不得宣称端到端成功。**

### 7.1 单元测试（先跑，最便宜）

```bash
cd /public/home/proj_kcchow/zhaoy/WRF_build/macda-to-wrfinterm
python -m unittest discover -s tests -v
```

检查点：`tests.test_mars_time` 与 `tests.test_converter` 全部通过；
不需要 `output/`、`public/` 或外部 NC 文件。

### 7.2 dry-run：先看标签，再写文件

```bash
python run_macda2w.py -c config/config.MACDA-v2-MY34-MY35.ini --dry-run
```

检查点：

- `label_anchor=time[0] +0034-12-53T02:00:00A`
- `records=360 relabelled=313 first=0034-00665_02:00:00 last=0035-00026_00:00:00`
- `times:` 表逐条打印，第 46/47/59/359 行与本文第 4 节表格一致
- 没有 `label spacing does not match` / `not strictly increasing` / `overwrite` 报错

### 7.3 锚点稳定性（部分重跑必须与全量同名）

```bash
python run_macda2w.py -c config/config.MACDA-v2-MY34-MY35.ini \
    --time-index 47 --dry-run
python run_macda2w.py -c config/config.MACDA-v2-MY34-MY35.ini \
    --time-index 0 --time-index 47 --dry-run
```

检查点：两条命令给索引 47 的标签都必须是 `0034-00669_00`，`xfcst=94 h`。
若出现 `0035-00001_00`，说明锚点又退化成了“当前选择的第一条”。

### 7.4 完整转换

```bash
python run_macda2w.py -c config/config.MACDA-v2-MY34-MY35.ini 2>&1 | tee c2w.MY34-MY35.log
ls output/MY34-MY35 | wc -l          # 期望 360
ls output/MY34-MY35 | sort | head -1 # 期望 MACDA:0034-00665_02
ls output/MY34-MY35 | sort | tail -1 # 期望 MACDA:0035-00026_00
ls output/MY34-MY35/MACDA:0034-00669_00   # 幻影 sol 必须存在
ls output/MY34-MY35/MACDA:0034-00669_22
```

检查点：360 个文件；`MACDA:0035-00027_00` **不应**存在；日志中每个文件都通过
写后回读校验（HDATE、维度、`XFCST`、必需字段齐全）。

> 注意：现有的 `output/MY34-MY35` 是早先混合口径的输出，**不能**当作验证结果。
> 请写入一个全新的空目录，或在确认无用后由使用者自行清理。

### 7.5 metgrid

`namelist.wps` 的 `&share` 必须与新标签一致：

```
 start_date       = '0034-00665_02:00:00',
 end_date         = '0035-00026_00:00:00',
 interval_seconds = 7200,
```

`&metgrid` 的 `fg_name` 需指向上一步的输出前缀（`MACDA`）。

```bash
cd $WPS
./metgrid.exe >& log.metgrid
grep -i "Couldn't open\|ERROR\|Aborting" log.metgrid
ls met_em.d01.* | wc -l              # 期望 360
```

检查点：

- **无** `Couldn't open file MACDA:0034-00669_00`（这正是本方案要消除的报错）
- `end_date` 若仍写 `0035-00027_00:00:00`，metgrid 会多要 2 个时次并报错；
  718 h / 2 h = 359 个间隔 → 360 个时次，末端标签就是 `0035-00026_00:00:00`
- `met_em.d01.0034-00669_00:00:00.nc` 存在，且 `met_em` 文件数为 360

### 7.6 real.exe

`namelist.input` 的 `&time_control` 起止时间与 `metgrid` 一致（起 MY 34 sol 665
02:00，止标签 MY 35 sol 26 00:00），`interval_seconds = 7200`。

```bash
cd $WRF/run
./real.exe >& log.real          # 或 mpirun -np N ./real.exe
grep -i "SUCCESS\|ERROR\|FATAL" log.real
ls -l wrfinput_d01 wrfbdy_d01
```

检查点：`SUCCESS COMPLETE REAL_EM INIT`；`wrfbdy_d01` 覆盖全部 359 个边界区间；
无日期越界或 `date not found` 类报错。

### 7.7 wrf.exe

```bash
cd $WRF/run
mpirun -np N ./wrf.exe >& log.wrf
grep -i "SUCCESS\|ERROR\|FATAL\|cfl" log.wrf | head
ncdump -v Times wrfout_d01_* | head -20
```

检查点：`SUCCESS COMPLETE WRF`；`Times` 序列连续、等间距；
跨越 `0034-00669` 处无中断。读取结果时按第 5 节把标签换回真实日期。

### 7.8 验证状态

| 项目 | 状态 |
| --- | --- |
| 单元测试 | **未执行**（本地禁止运行） |
| dry-run | **未执行** |
| 完整转换 360 文件 | **未执行** |
| metgrid | **未执行** |
| real.exe | **未执行** |
| wrf.exe | **未执行** |

上表全部转为“通过”之前，本方案只能被描述为“已完成代码实现、待 HPC 验证”，
**不得**声称跨年模拟已端到端成功。
