"""MCD background entry point. Use --help; no log files are generated."""

# 命令行入口只组织处理步骤，不在这里重复实现时间、坐标或沙尘公式。
# 顺序为解析参数、读取配置和源坐标、选择时段、写全年日历、按 sol 写背景，最后可选检查/导出。
# 默认生成月平均插值背景；--eof 生成合成扰动候选，两种输出都没有观测同化。

# argparse 处理命令行，configparser 读取 ini，groupby 将连续时刻归入同一个 sol。
import argparse
import configparser
from dataclasses import asdict
from itertools import groupby
import json
from pathlib import Path

# 复用物理参数对象和各模块接口，使命令行与测试调用同一套计算函数。
from macda2wrf.mcd_dust import DustBinParameters
from macda2wrf.mcd_pipeline import write_calendar, write_states, export_intermediate
from macda2wrf.mcd_reader import McdReader
from macda2wrf.mcd_time import iter_mcd_times
from macda2wrf.mcd_value import McdVariability


# 配置和变量表相对脚本所在项目定位，运行目录改变也不会误用另一个同名文件。
PROJECT_DIR = Path(__file__).resolve().parent


def arguments():
    # description 使用原模块 docstring；新增中文注释不会改变 --help 的现有输出。
    parser = argparse.ArgumentParser(description=__doc__)
    # 默认 MY36 ini 可被 --config 替换；start/end-sol 在年内从 1 开始编号，端点均包含。
    parser.add_argument("--config", type=Path, default=PROJECT_DIR / "conf/config.MCD-MY36.ini")
    parser.add_argument("--start-sol", type=int)
    parser.add_argument("--end-sol", type=int)
    # max-times 只缩短所选物理场的处理列表，不缩短独立日历的全年覆盖。
    parser.add_argument("--max-times", type=int, help="Explicitly limit processing; calendar still spans full MY")
    # calendar-only 只生成/核验日历；dry-run 只读配置和必要坐标、打印计划，不写产品。
    parser.add_argument("--calendar-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    # EOF 扰动须显式启用，并写入独立候选目录；不能把它理解为每日真值恢复开关。
    parser.add_argument("--eof", action="store_true", help="Optional synthetic coherent PCM EOF variability, not truth")
    # WPS 输出也须显式启用；missing 保留超范围缺测并使导出拒绝，endpoint 允许端点延拓。
    # 两种选项都不允许向当地地表以下外推；尘埃 TRC 字段还有单独的守恒接收限制。
    parser.add_argument("--write-intermediate", action="store_true")
    parser.add_argument("--wps-extrapolation", choices=("missing", "endpoint"), default="missing")
    # validate 做独立物理检查；resume 只复用已通过检查且配置/日期匹配的已有 sol 文件。
    parser.add_argument("--validate", action="store_true", help="Independently check physical NetCDF output")
    parser.add_argument("--resume", action="store_true", help="Check and reuse matching complete sol files")
    return parser.parse_args()


def main():
    # 先解析命令行，再读取 UTF-8 配置；配置文件不存在就立即停止，避免沿用空配置。
    args = arguments()
    config = configparser.ConfigParser()
    if not config.read(args.config, encoding="utf-8"):
        raise FileNotFoundError(args.config)
    # 初始化读取器只建立月资料网格、混合层和地形基准；完整月场在实际采样时再加载。
    source_dir = Path(config["paths"]["source_dir"])
    output_dir = Path(config["paths"]["output_dir"])
    reader = McdReader(source_dir, config["paths"]["mountain_file"], config.getint("time", "year"))
    # ini 中 dust 键与 DustBinParameters 字段一一对应；转浮点后由参数类核验合法范围。
    parameters = DustBinParameters(**{name: float(value) for name, value in config["dust"].items()})
    step_hours = config.getint("time", "step_hours")
    if step_hours != 2:
        raise ValueError("this MY36 product requires a two-Mars-hour interval")
    # 时间模块统一平/真太阳时、Ls、SI 秒和 MarsWRF 日期；2h 指两个火星标签小时。
    # MY36 标签年包含 669 sols，全年为 8028 个时刻；668.6-sol 理想轨道末尾另带回绕标志。
    all_points = list(iter_mcd_times(reader.year, step_hours))
    # 命令行覆盖 ini 的 sol 范围；未指定的项沿用配置，不修改原配置文件。
    start_sol = args.start_sol if args.start_sol is not None else config.getint("time", "start_sol")
    end_sol = args.end_sol if args.end_sol is not None else config.getint("time", "end_sol")
    if not 1 <= start_sol <= end_sol <= all_points[-1].sol_of_year:
        raise ValueError("requested sol interval is outside the model year")
    # 先按闭区间筛选 sol，再按 max-times 截断，适合小样本检查或明确限制处理范围。
    points = [point for point in all_points if start_sol <= point.sol_of_year <= end_sol]
    if args.max_times is not None:
        if args.max_times < 1:
            raise ValueError("max-times must be positive")
        points = points[:args.max_times]
    # 标准输出给出所选/全年样本数、日期、网格与目录；此脚本不创建日志文件。
    print(json.dumps({"selected_times": len(points), "full_year_times": len(all_points),
        "first": points[0].wrf_date, "last": points[-1].wrf_date,
        "grid": [len(reader.latitude), len(reader.longitude)], "source": str(source_dir),
        "output": str(output_dir), "synthetic_eof": args.eof}, ensure_ascii=False))
    if args.dry_run:
        return
    # 日历始终使用全年 all_points；因此截取几天处理也仍能保留统一完整的时间基准。
    calendar_path = output_dir / f"MY{reader.year}_calendar_2h.nc"
    if not calendar_path.exists():
        write_calendar(calendar_path, reader, all_points, parameters)
        print(f"calendar: {calendar_path}")
    else:
        # 已有日历不覆盖，只比对长度和全部日期标签；这里未声称核验其所有属性。
        from netCDF4 import Dataset

        with Dataset(calendar_path) as calendar:
            if len(calendar.dimensions["time"]) != len(all_points):
                raise ValueError("existing calendar is incompatible; no file overwritten")
            if list(calendar["wrf_date"][:]) != [point.wrf_date for point in all_points]:
                raise ValueError("existing calendar date strings are incompatible")
    if args.calendar_only:
        return
    variability = None
    if args.eof:
        # seed 控制可复现的扰动序列，amplitude 控制幅度，modes 限制使用的 EOF 模态数。
        variability = McdVariability(config["paths"]["eof_file"],
            seed=config.getint("variability", "seed"), amplitude=config.getfloat("variability", "amplitude"),
            modes=config.getint("variability", "modes"))
    # 用独立目录区分合成候选和未扰动背景，避免下游诊断把二者混为同一产品。
    kind = "eof_candidate" if args.eof else "background"
    destination = output_dir / kind
    variable_table = PROJECT_DIR / config["paths"]["variable_table"]
    # WPS 目标压力层由逗号分隔，单位 Pa；真正使用时还会检查单调性和地下层。
    pressure_levels = [float(value) for value in config["wps"]["pressure_levels_pa"].split(",")]
    for sol, group in groupby(points, key=lambda point: point.sol_of_year):
        # points 已按时间排序，groupby 每次拿出同一个 sol；文件名含起止小时，表示本块覆盖。
        block = list(group)
        tag = f"MY{reader.year}_sol{sol:05d}_{block[0].date.hour:02d}-{block[-1].date.hour:02d}_{kind}"
        path = destination / f"{tag}.nc"
        if path.exists() and args.resume:
            # 续算不是只看文件是否存在：先做物理校验，再核对影响结果的参数、源目录和时序。
            from netCDF4 import Dataset
            from utils.validate_mcd import validate

            report = validate(path)
            if not report["valid"]:
                raise ValueError(f"existing file failed validation: {path}")
            with Dataset(path) as existing:
                # JSON 排序后逐字比较参数，避免续算时偷偷混用另一套粒径或光学假设。
                if existing.dust_parameters != json.dumps(asdict(parameters), sort_keys=True):
                    raise ValueError(f"existing dust parameters differ: {path}")
                if existing.source_directory != str(reader.source_dir):
                    raise ValueError(f"existing source directory differs: {path}")
                if list(existing["wrf_date"][:]) != [point.wrf_date for point in block]:
                    raise ValueError(f"existing time coverage differs: {path}")
                if args.eof:
                    # 合成候选还必须匹配随机种子、幅度、请求模态数和 EOF 文件路径。
                    metadata = json.loads(existing.variability_metadata)
                    requested = {"seed": variability.seed, "amplitude": variability.amplitude,
                        "requested_modes": variability.requested_modes,
                        "source_file": str(variability.eof_path.resolve())}
                    if any(metadata.get(name) != value for name, value in requested.items()):
                        raise ValueError(f"existing EOF settings differ: {path}")
            print(f"checked existing state: {path}", flush=True)
        else:
            # 未启用 resume 时，即使路径存在也交给 write_states 的拒绝覆盖检查处理。
            write_states(path, reader, block, parameters, variability)
        if args.validate:
            # 独立检查可以作用于新文件或续算复用文件；未通过即停止后续 sol 和 WPS 导出。
            from utils.validate_mcd import validate

            report = validate(path)
            print(json.dumps(report, ensure_ascii=False))
            if not report["valid"]:
                raise ValueError(f"validation failed: {path}")
        if args.write_intermediate:
            for point in block:
                # 每个时刻单独导出 WPS 文件，文件名复用火星 YYYY-DDDDD_HH:MM:SS 标签。
                # 这里仍会重新构建该时刻场，且不宣称 WPS、real.exe 或 WRF 已成功读取。
                target = output_dir / "intermediate" / f"MCD:{point.wrf_date}"
                export_intermediate(target, reader, point, pressure_levels, variable_table,
                    parameters, extrapolation=args.wps_extrapolation, variability=variability)
        print(f"state: {path}", flush=True)


if __name__ == "__main__":
    # 直接执行脚本才启动转换；被测试或其他模块 import 时不触发写文件。
    main()
