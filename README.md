## 问题探索

MarsWRF 普通 Mars 模式并不在 I/O 时把地球日期动态转换成火星日期，而是通过 `WRF_PLANET` 编译为独立日历：

```text
YYYY-DDDDD_HH:MM:SS
```

实际文件验证结果：

- `wrfinput_d01`: `0001-00001_00:00:00`
- `wrfout`: `0003-00001_00:00:00`
- `wrfrst`: `0001-00520_00:00:00`
- `PLANET_YEAR=669`
- `P2SI=1.027491`，一个物理 sol 约 88775.2 SI 秒

只有 `MARS24_TIMING` 模式才通过 UTC Julian day、TT、J2000、MSD 进行天文时间转换。

完整源码说明和 `/phys` 火星代码块清单见 [MARSWRF_AUDIT.md](/home/zy/WRF/record-for-atmos-geos/macda-to-wrfinterm/MARSWRF_AUDIT.md:1)。

## 时间修复

新增了 [mars_time.py](/home/zy/WRF/record-for-atmos-geos/macda-to-wrfinterm/macda2wrf/mars_time.py:20)，实现 MACDA 的 `669/668/669/668/669` 五年周期和月份长度。

目标文件首时次现在正确转换为：

```text
+0028-10-07T02:00:00A
time=3180.08333333333 sol
-> MY28 SOY507
-> 0028-00507_02:00:00.0000
```

连续 `time` 和 `Mars_date` 会相互校验，不一致时停止转换。

## 变量修复

[MACDA-v2_CORE.csv](/home/zy/WRF/record-for-atmos-geos/macda-to-wrfinterm/db/MACDA-v2_CORE.csv:1) 已修正：

- `temp/uwind/vwind/geop` 转压力层 `TT/UU/VV/GHT`
- 气压源字段改为 `PRESSURE`，由 metgrid 派生最终 `PRES`
- `PSFC`、`SKINTEMP` 直接输出
- `SPECHUMD/QV` 暂置零
- `TAU_OD2D` 归一到 700 Pa
- `CO2ICE` 直接输出
- `UU/VV` 正确标记为地理坐标相对风

`omega` 不能直接等同 WRF 几何垂直速度；`dustmmr` 也不能在没有粒径分配假设时直接拆成 `TRC01/TRC02`，因此没有做错误映射。

## 验证

真实首、次、末时次均已生成，每个文件包含 137 条记录、72×36 网格；末时次 `XFCST=718` Martian hours。7 个单元测试通过，WPS 4.6 `rd_intermediate.exe` 也成功读取全部记录。

需要注意：旧的 `marswrf/WPS` 日期工具尚未包含五位 sol 日历分支，而且缺少 Mars 静态地理数据，所以当前已完成并验证的是 `MACDA -> WRF intermediate`；`metgrid -> real.exe` 还不能据此宣称完整可运行。默认 Python 环境还需按 [requirements.macda.txt](/home/zy/WRF/record-for-atmos-geos/macda-to-wrfinterm/requirements.macda.txt:1) 安装 `h5py`。