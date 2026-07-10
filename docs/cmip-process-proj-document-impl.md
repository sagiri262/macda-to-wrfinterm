# cmip6-to-wrfinterm
### *.ini 的作用
conf/config.<MODEL>.ini 配置输入目录、模型名、scenario/experiment、ensemble/grid 标识、ETL 起止时间、输出目录等。入口脚本按 -m 找配置文件。
output_prefix 这类字段现在基本不被主流进程所采用，输出的前缀来自于 [INPUT] model_name。
### /db/*.csv 文件
csv文件主要内容的功能类似于 Vtable，包括
```text
src_v, aim_v, units, type, lvlmark, desc
```
也就是源变量名、WRF intermediate 字段名、单位、变量类型、层次标记和描述。其中 type 很关键，它决定代码走 3D、2D、土壤、月/日频率、固定场等哪条处理逻辑。

还要补充一个核心表：db/cmip6_meta.csv。它决定某个模型和 scenario 要加载哪些变量组、对应哪个 CMIP6 table、主输出频率是哪一个。

### /lib 目录

**lib/cmip_handler.py** 是核心调度层：读 meta/vtable、加载 NetCDF、选时间、插值、缺测填补、垂直层处理、土壤处理、写 intermediate 文件。
**lib/adapters/*.py** 负责模型相关的文件发现和打开方式。标准 CMIP6 文件名模式是：
```text
<var>_<table>_<model>_<scenario>_<member>_<grid>[_<time-range>].nc
```
**utils/grid.py** 定义输出 intermediate 的规则经纬网、气压层、土壤层，默认是全球 1 度、181x360、14 个气压层。
它不是读取 WRF namelist.input 里的 dx/dy 来重采样到 WRF 网格。当前代码是把 CMIP 数据插值到 OutputGrid 的规则经纬度网格，默认全球 1 度；如果要改，需要在 .ini 的 [OUTPUT] 里配置 lat_start/lon_start/nlat/nlon/plev_hpa 等。最后的插值到WRF投影网格，是后续 ./metgrid.exe 程序执行的，包括插值到 2D 平面网格和 Lev 垂直网格。

**单位换算**是通过 .**/lib/cmip_handler.py** 文件实现。单位换算不是通用单位系统，代码里有少量硬编码单位转换，比如 tos 从摄氏度转 K、mrsos 的土壤湿度处理。它不是一个完整的“自动单位换算框架”。

### sample 目录
sample 不是只放测试输入数据。它还包含下载脚本、生成假数据脚本、namelist.wps、namelist.input。比如 CESM2 SSP245 样例说明了下载、转换、链接到 WPS 的流程。

### run_c2w.py 
该程序是主入口。它读取配置，构造 CMIPHandler，遍历输出时间序列，调用 parse_data() 和 write_wrfinterm()。它不是“调用/模拟 ungrib.exe 的输出过程”，而是直接写 WPS intermediate format 文件。效果上，这些文件扮演了 ungrib.exe 输出的角色，后续可以直接给 metgrid.exe 读。也就是说通常流程是：

```text
CMIP6 NetCDF -> run_c2w.py -> WRF intermediate files -> metgrid.exe -> met_em -> real.exe/wrf.exe
```

## 总结
这个项目用 .ini 选择数据源和时间段，用 cmip6_meta.csv 选择变量组和频率，用各模型 CSV/Vtable 定义变量语义，用 adapter 找 NetCDF 文件，用 CMIPHandler 做时间选择、水平插值、垂直层/土壤处理和少量单位处理，最后直接写出可被 metgrid.exe 读取的 WRF intermediate 文件。
