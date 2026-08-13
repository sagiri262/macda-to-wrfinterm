# MACDA 沙尘初始化

## 转换约定

`db/MACDA-v2_CORE.csv` 将必需的三维 MACDA `dustmmr` 映射为 MarsWRF 标量
`TRC01` 和 `TRC02`：

```text
TRC01 = dustmmr * dust_bin1_mass_fraction / marswrf_wt1
TRC02 = dustmmr * dust_bin2_mass_fraction / marswrf_wt2
```

默认质量分数和辐射权重均为 `0.3/0.7`。因此两个存储 tracer 都等于 `dustmmr`，但辐射使用的
总质量严格守恒：

```text
marswrf_wt1 * TRC01 + marswrf_wt2 * TRC02 = dustmmr
```

转换器会在每个配置压力层写出两个字段，并回读输出。
缺少压力层、出现负值或加权沙尘全零时，转换会失败。
加权结果无法重建插值后的源 `dustmmr` 时，转换同样会失败。
源文件缺少 `dustmmr` 也会报错，不能静默产生无沙尘强迫文件。

## HPC 接力检查

将 `sample/metgrid/METGRID.TBL.MARS.dust` 中 `TRC01` 和 `TRC02` 的条目加入
`opt_metgrid_tbl_path` 实际指向的 `METGRID.TBL`。只修改未启用副本不会生效。

在独立 HPC 案例目录中按顺序执行：

```bash
python run_macda2w.py --config config/config.MACDA-v2.ini
./metgrid.exe
ncdump -h met_em.d01.<time>.nc | grep -E 'TRC01|TRC02'
ncdump -v TRC01,TRC02 met_em.d01.<time>.nc
./real.exe
ncdump -h wrfinput_d01 | grep -E 'TRC01|TRC02'
ncdump -v TRC01,TRC02 wrfinput_d01
```

| 检查点 | 正常结果 | 缺失或全零时的位置 |
| --- | --- | --- |
| Intermediate `MACDA:*` | 两个 tracer 存在且至少部分网格非零 | 转换器或使用了旧 intermediate |
| `met_em` | 两个 tracer 存在且至少部分网格非零 | 实际 `METGRID.TBL` 没有传递字段 |
| `wrfinput_d01` | 两个 tracer 存在且至少部分网格非零 | `real.exe` 读入或初始化环节 |
| 首个 `wrfout` 时次 | 物理积分前两个 tracer 非零 | 初始化后、输出前被覆盖 |

`mp_physics = 48` 的 WRF Registry 将两个沙尘标量状态变量命名为 `TRC01` 和 `TRC02`。
现有 `module_initialize_real.F` 审计未发现它们的专用压力层输入分支。因此 `met_em` 非零而
`wrfinput_d01` 为零时，问题在 WRF 端，不是转换器或地表库存。

## 地表沙尘开关

`reset_surfdust1/2` 在模式启动时重置 `surftracer1/2`，影响之后在线起沙是否有地表库存；
它们不读入、不写出、也不初始化大气 `TRC01/TRC02`。验证 MACDA 初始化时，应先通过上表四个
检查点确认大气 tracer 非零，再单独进行在线起沙和沉降实验。

## 光学厚度

MACDA v2.0 的 `dustmmr` 使用有效半径 1.5 微米、颗粒密度 2500 kg m-3 和 637 nm
`Qext = 2.62718`。当前 MarsWRF 配置为 `ddpp = 1.0e-6 m`、`ddpp2 = 2.4e-6 m`、
颗粒密度 2650 kg m-3、`Qext1 = Qext2 = 3.04` 和权重 0.3/0.7。质量守恒时，预计
MarsWRF/MACDA 光学厚度比例约为 1.934。这是粒径档光学参数不匹配，不表示大气沙尘缺失。

本地输入文件的首时次中，按 MACDA 参数由 `dustmmr` 重建的光学厚度与 `coldust` 的中位比为
0.991825；按当前 MarsWRF 双粒径参数则为 `coldust` 的 1.9182 倍。若需光学厚度闭合，应单独
标定双粒径档参数，不能破坏上面的质量守恒关系。
