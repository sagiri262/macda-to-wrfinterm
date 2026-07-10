## 对齐变量名

查询 mro-mcs* 文件的头文件等信息，
HDF5 "mro-mcs-reanalysis_mars_MY28SOY237_MY28SOY267_v2-0.nc" {
FILE_CONTENTS {
 group      /
 dataset    /Ls
 dataset    /MY_Ls
 dataset    /Mars_date
 dataset    /co2ice
 dataset    /coldust
 dataset    /dustmmr
 dataset    /geop
 dataset    /lat
 dataset    /lev
 dataset    /lon
 dataset    /lwflux
 dataset    /omega
 dataset    /psurf
 dataset    /swflux
 dataset    /temp
 dataset    /time
 dataset    /tsurf
 dataset    /uwind
 dataset    /vwind
 }
}

**一一确认每一个变量对应的内容**
Ls ——> 火星太阳黄经
MY_Ls ——> 基于火星太阳黄经划分的火星年份编号，用于跨年比较同一季节。
Mars_data ——> 火星本初子午线处的日期/时间
co2ice ——> 地表二氧化碳冰柱质量
coldust ——> 柱积分尘埃光学深度
dustmmr ——> 尘埃质量混合比
geop ——> 位势
lat ——> 纬度
lon ——> 精度
lev ——> 垂直 sigma 模式层，无量纲地形跟随坐标。到 WRF 需要转换到 eta 坐标
lwflux ——> 到达地表的向下长波辐射通量
omega ——> 压力坐标中的垂直速度，即气压随气团运动的变化率。小于零表示下沉，大于零表示上升
psurf ——> 地表气压
swflux ——> 到达地表的向下长波辐射通量
temp ——> 大气温度
time ——> 数据时间坐标。CEDA 标为 days，同时数据组织上每文件覆盖 30 sols、每 sol 12 时次.
tsurf ——> 地表温度
uwind ——> u 分向风速
vwind ——> v 分向风速


下面列出完整的，WRF 所需要的最小输入变量的名称、物理意义和变量描述。但是要特别注意，比如 Relative Humidity 在火星上是没有的，可能是需要单独计算或者干脆不输入的；再比如土壤温湿度更不可能有了，MarsWRF中也说了，这个单独置0，不输入到WRF里面去。
![alt text](image.png)


```
temp,TT
uwind,UU
vwind,VV
geop,GHT
psurf,PSFC
tsurf,SKINTEMP
```

整理可以直接输入到 MarsWRF 里的变量大概有这些，对照 CESM2_3hr.csv 的内容，我们大概可以写成这个形式

src_v,aim_v,units,type,lvlmark,desc
temp,TT,K,2d,PlevPt,2-m temperature
uwind,UU,m s-1,2d,PlevPt,10m wind u-component
vwind,VV,m s-1,2d,PlevPt,10m wind v-component
psurf,PSFC,Pa,

## 修改程序
参考项目 cmip6-to-wrfinterm 。主要需要修改的代码