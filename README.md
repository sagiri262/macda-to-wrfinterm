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

# MarsWRF 的基本流程

整体流程是：
```text
MACDA NetCDF
    |
    | run_macda2w.py
    v
MACDA:0028-00507_02 等 WRF intermediate 文件
    |
    |                       火星静态地理 tile
    |                              |
    |                        geogrid.exe
    |                              |
    |                        geo_em.d01.nc
    |                              |
    +----------- metgrid.exe <-----+
                    |
                    v
       met_em.d01.0028-00507_02:00:00.nc
                    |
                 real.exe
                    |
             wrfinput / wrfbdy
                    |
                 wrf.exe
```

首先，地球WPS不能直接跑，要先修改WPS里的内容。在 /docs 目录下的 `验证可行性指南` 里面已经讲清楚了，哪些东西是我们需要改的。改完之后，就要重新编译一下WPS里的可执行程序，主要是 ./geogrid.exe 和 ./metgrid.exe。
由于修改了半径常数，需要在 HPC 环境中重新编译：
```shell
cd /HPC路径/marswrf/WPS

./clean -a
./configure

# 其实只需要执行这两个就可以了
# 如果不行再重新编译整个WPS
./compile geogrid
./compile metgrid
```

# 生成 MACDA intermediate

安装依赖：
```shell
cd /HPC路径/marswrf/macda-to-wrfinterm

python -m pip install -r requirements.macda.txt
```

先检查：
```shell
python run_macda2w.py \
  -c config/config.MACDA-v2.ini \
  --dry-run
```

生成全部 360 个时次：
```shell
python run_macda2w.py \
  -c config/config.MACDA-v2.ini \
  --max-times 360
```

不能只执行不带参数的默认命令，因为当前配置中 `max_times = 1` 默认只会生成第一个时次。
输出类似：
```text
output/MACDA:0028-00507_02
output/MACDA:0028-00507_04
...
output/MACDA:0028-00537_00
```

这些是 WRF intermediate，不是 geo_em、met_em 或 wrfinput。当前不建议用 HPC job array 分别执行单个 --time-index，因为每个独立进程都会把自己的首时次作为 XFCST=0，导致整段序列的 XFCST 不连续。先用一次命令顺序生成全部时次。

# 执行 geogrid
在独立 WPS 运行目录中准备：
```textnamelist.wps
geogrid.exe
geogrid/GEOGRID.TBL
WPS_GEOG_MARS
```

然后执行：
```shell
./geogrid.exe
```

预期生成`geo_em.d01.nc`文件，geogrid 只需要在网格范围、投影、分辨率或静态地理数据改变时重新执行，不需要每个 MACDA 时次执行一次。

# 执行 metgrid
将 intermediate 文件放在 WPS 运行目录，或者建立链接 `ln -s /HPC路径/macda-to-wrfinterm/output/MACDA:* .`，并且要求确保 `namelist.wps` 里的这些内容修改到位：
```fortran
&ungrib
 prefix = 'MACDA'
/

&metgrid
 fg_name = 'MACDA'
/
```

还要准备与字段映射一致的 `metgrid/METGRID.TBL`，其重点是让 metgrid 从 PRESSURE 派生最终 PRES，接收 `TT/UU/VV/GHT/PSFC/SKINTEMP` 等WRF的输入变量，最后根据需要传递 TAU_OD2D 和 CO2ICE，但是这两个变量其实可有可无。
确保上面做的都没有什么问题后，执行一次 `./metgrid.exe`，metgrid 会自动遍历 start_date 到 end_date，按 interval_seconds=7200 处理全部时次，不需要手工对每个 intermediate 文件执行一次。
