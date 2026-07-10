# 0002 MarsWRF 时间与字段对齐

## 修改

- 增加有文献定义的 MACDA `669/668` sol 日历转换，并用连续 `time` 坐标做交叉校验。
- WRF intermediate 中的气压输入字段由 `PRES` 改为 `PRESSURE`；最终 `PRES` 由 metgrid 派生。
- 将 `UU` 和 `VV` 标记为相对地理坐标的风。
- 增加 NetCDF packed/missing value 处理，并加强配置文件和变量表校验。
- 增加 WRF intermediate 回读校验和单元测试。
- 替换空的样例 CSV 占位文件以及无效的样例 namelist 内容。
- 删除未被使用且只有重复表头的 `db/MACDA_LEV.csv`；`db/MACDA-v2_CORE.csv` 成为唯一生效的 MACDA 变量表。
- 记录 MarsWRF 时间实现、实际 WRF/MACDA 元数据、`phys` 中火星相关代码块以及下游 WPS 限制。

复制项目中实际去掉的旧占位 CSV 为：

```text
db/MACDA_LEV.csv
sample/MACDA-v2/db/MACDA_ATM.csv
sample/MACDA-v2/db/MACDA_SFC.csv
sample/MACDA-v2/db/macda_meta.csv
```

这些文件没有被 `run_macda2w.py` 或 `macda2wrf` 包读取，而且内容为空或只是重复表头。继续保留会让人误以为存在多套生效映射，因此统一改由 `db/MACDA-v2_CORE.csv` 管理字段。

## 已核对的转换

```text
+0028-10-07T02:00:00A
time=3180.08333333333 sol
-> 0028-00507_02:00:00.0000
```

