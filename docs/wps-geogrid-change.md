| geogrid 字段 | 推荐来源 | 处理方式 |
|---|---|---|
| `HGT_M` | MGS/MOLA MEGDR 地形 | 转成 WPS geogrid binary tile |
| `THC/THCBCK` | MGS/TES thermal inertia 或相关论文产品 | 插值到 geogrid 静态场 |
| `ALBEDO/ALBBCK` | TES/OMEGA/MCD 类全球反照率产品 | 先常数，后替换真实图 |
| `EMISS/EMBCK` | TES emissivity 或常数 | 初期可设 0.95/1.0 |
| `H2OICE/GRD_ICE_*` | Mars ice/public climatology | 后续增强项 |
| dust | MACDA `coldust` 或 Montabone dust climatology | 更适合做时变 intermediate，不是 geogrid 静态场 |

公开资料入口：WPS 官方说明见 https://www2.mmm.ucar.edu/wrf/users/wrf_users_guide/build/html/wps.html ，MOLA MEGDR 见 https://pds-geosciences.wustl.edu/missions/mgs/megdr.html ，MGS/TES 资料见 https://pds-geosciences.wustl.edu/missions/mgs/tes.html ，TES dust/ice 光学厚度说明见 https://atmos.nmsu.edu/PDS/data/PDS4/mgs_tes_atmos_dust-ice/document/user_guide_TES_COD.pdf 。