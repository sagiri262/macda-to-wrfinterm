# 0003 WPS geogrid 文件准备

## 新增

- `sample/namelist.wps.copy`
- `sample/geogrid/GEOGRID.TBL.ARW.copy`
- `sample/geogrid/GEOGRID.TBL.copy`
- `geogrid/GEOGRID.TBL.ARW.copy`
- `sample/geogrid/GEOGRID.TBL.WPS-4.6.0.reference.copy`
- `sample/geogrid/geogrid.exe.copy`
- `sample/geogrid/lib/libnetcdff.so.7.1.0.copy`
- `sample/geogrid/lib/libnetcdf.so.19.copy`
- `sample/WPS_GEOG_MARS/*/index.copy` 共 11 份
- `docs/GEOGRID_FILE_AUDIT.md`
- `docs/GEOGRID_COPY_CHANGES.md`
- `docs/GEOGRID_MANUAL_VALIDATION.md`

## 修改

没有修改任何原文件。所有配置结果都使用 `.copy` 后缀。

## 设计结论

- 网格对齐为单域全球 5 度 `lat-lon`，WPS 定义 `73 x 37`，对应 MACDA `72 x 36` mass grid。
- 时间对齐为 MY28/SOY507 02:00 至 SOY537 00:00，间隔 7200 个火星时钟秒。
- geogrid 表移除地球植被、雪、城市等静态场，保留火星地形、全陆面分类、反照率、热惯量、发射率和可选水冰。
- 已有 MOLA/TES 原始文件不是 WPS tile，因此只建立 index 模板，没有复制大文件或伪造数据格式。
- 现成 WPS 4.6 二进制硬编码地球半径，旧 Mars WPS 日期代码不支持五位 sol；二者均作为后续源码修改和重新编译前的明确限制保留。

## 执行范围

按任务要求，没有编译或运行 geogrid、metgrid、real 或 wrf。
