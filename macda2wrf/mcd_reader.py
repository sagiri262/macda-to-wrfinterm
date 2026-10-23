"""Read and interpolate native MCD monthly NetCDF3 fields, without regridding."""

# 本模块负责把月文件读成“指定季节、指定本初子午线真太阳时”的全球背景场。
# 处理顺序是读取源坐标与地形、缓存相邻月份、做日变化插值，再做季节插值。
# 水平网格和混合层索引均保持源资料定义；这里不做 WRF 网格或 eta 层转换。
# 插值给出连续的气候背景，不能从月平均资料恢复某个真实 sol 的天气。


from collections import OrderedDict
from pathlib import Path
from netCDF4 import Dataset
import numpy as np
from .mcd_time import diurnal_weights


# 月平均文件的源变量名；三维大气量带 Time/层/纬度/经度，地表量带 Time/纬度/经度。
SOURCE_FIELDS = ("temp", "u", "v", "rho", "dustq", "reffdust",
    "vmr_h2ovapor", "ps", "tsurf", "co2ice", "h2oice", "tau_pref_gcm")
# sd 文件提供无符号 RMS 变化尺度；改名时保留 rms_ 前缀，防止被当作有符号距平。
RMS_FIELDS = {"rmstemp": "rms_temperature", "rmsu": "rms_u", "rmsv": "rms_v",
    "rmsps": "rms_ps", "rmstau_pref_gcm": "rms_tau_visible610_mcd"}
# 这些字符串是源文件单位契约，用于核验而不是自动换算；拼写不符就拒绝读取。
# dustq 已经是沙尘质量/空气质量，reffdust 已经以米表示，NU 表示源光学厚度无量纲。
EXPECTED_UNITS = {"temp": "K",
                  "u": "m.s-1",
                  "v": "m.s-1",
                  "rho": "kg/m3",
                  "dustq": "kg/kg",
                  "reffdust": "m",
                  "ps": "Pa",
                  "tsurf": "K",
                  "co2ice": "kg.m-2",
                  "h2oice": "kg/m2",
                  "vmr_h2ovapor": "mol/mol",
                  "tau_pref_gcm": "NU"}


def read_array(variable):
    """Decode NetCDF masks into explicit NaNs instead of exposing fill values."""
    # 先把掩码对应的缺测值换成 NaN，再统一为双精度；不让 _FillValue 混进物理运算。
    return np.asarray(np.ma.filled(variable[:], np.nan), dtype=np.float64)


class McdReader:
    """
    Cache at most two complete months; each sample is a simultaneous map.

    MCD Time is prime-meridian true solar time. Do not add longitude to this
    index: its spatial maps already describe one simultaneous global state.
    The year boundary policy deliberately reuses the selected year's months.
    It is a labelled climatological cycle, not adjacent-year weather.
    """

    def __init__(self, source_dir, mountain_file, year=36):
        # 保存月资料目录、独立地形文件和火星年；此处不创建或修改任何源文件。
        self.source_dir = Path(source_dir)
        self.mountain_file = Path(mountain_file)
        self.year = year
        # 缓存值是已经读入内存的数组字典，不是需要长期保持打开的 NetCDF 句柄。
        self.cache = OrderedDict()
        # 用第一个月建立坐标基准。with 结束即关闭文件，后续月份必须符合这份基准。
        with Dataset(self.month_path(1)) as source:
            self.latitude = read_array(source["latitude"])
            self.longitude = read_array(source["longitude"])
            # 层中心压力由 p(k,lat,lon)=aps(k)+bps(k)*ps(lat,lon) 重建。
            self.aps = read_array(source["aps"])
            self.bps = read_array(source["bps"])
            # altitude 是源文件的高度标签，保存用于溯源；不把它直接当几何高度或 GHT。
            self.source_altitude_km = read_array(source["altitude"])
            # Time 为本初子午线真太阳时采样，通常为 2、4、…、24 火星时。
            # 每张图已代表同一物理时刻的全球场，不能对每个经度再次平移这个时间索引。
            self.hours = read_array(source["Time"])
        # 检查一维水平坐标是否等间距；这里允许源纬度为降序，输出方向由后续模块整理。
        for values in (self.latitude, self.longitude):
            if not np.allclose(np.diff(values), np.diff(values)[0]):
                raise ValueError("MCD horizontal grid is not regular")
        with Dataset(self.mountain_file) as source:
            self.validate_grid(source)
            # orography 是相对当地 areoid 的地形高；areoid 是行星基准面的半径，二者均为米。
            # 这两个量会用于静力高度重建，不能把 areoid 半径再次当作地形高度相加。
            self.orography = read_array(source["orography"])
            self.areoid_radius = read_array(source["areoid"])

    def month_path(self, month, suffix="me"):
        # 默认定位月平均 me 文件；调用者传入 sd 时定位同年、同月的 RMS 文件。
        return self.source_dir / f"MY{self.year}_{month:02d}_{suffix}.nc"

    def validate_grid(self, source):
        # 逐坐标比较长度与数值，拒绝把不同水平网格的数据按相同数组下标直接混合。
        # rtol=0 表示只允许 2e-5 度绝对舍入差，不随经纬度绝对值放宽容差。
        for name, expected in (("latitude", self.latitude), ("longitude", self.longitude)):
            actual = read_array(source[name])
            if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0, atol=2.0e-5):
                raise ValueError(f"{source.filepath()}: {name} does not match the native MCD grid")

    def load_month(self, month):
        # 命中缓存就把该月移到队尾，标记为最近使用；返回现有数组，避免重复磁盘读取。
        if month in self.cache:
            self.cache.move_to_end(month)
            return self.cache[month]
        fields = {}
        with Dataset(self.month_path(month)) as source:
            self.validate_grid(source)
            # 混合坐标系数和日内采样必须逐值相同，否则相同层/时刻下标没有共同含义。
            for name, expected in (("aps", self.aps), ("bps", self.bps), ("Time", self.hours)):
                if not np.array_equal(read_array(source[name]), expected):
                    raise ValueError(f"month {month}: coordinate {name} changes")
            for name in SOURCE_FIELDS:
                variable = source[name]
                # 先检查单位再读数值；未声明的单位差异不能靠变量名称猜测或静默修正。
                if variable.units != EXPECTED_UNITS[name]:
                    raise ValueError(f"{name}: unexpected source units {variable.units!r}")
                fields[name] = read_array(variable)
        # RMS 文件单独读取；当前代码只检查其水平网格和非负性，不额外核验其单位属性。
        with Dataset(self.month_path(month, "sd")) as source:
            self.validate_grid(source)
            for name, output_name in RMS_FIELDS.items():
                fields[output_name] = read_array(source[name])
                # RMS 是变化幅度，负值没有该定义下的意义；NaN 的最终拒绝由 build_state 完成。
                if np.any(fields[output_name] < 0):
                    raise ValueError(f"{name}: unsigned RMS must be nonnegative")
        self.cache[month] = fields
        # 相邻月插值最多需要两个月；移除最久未使用月份，限制全年处理的内存占用。
        while len(self.cache) > 2:
            self.cache.popitem(last=False)
        return fields

    def sample(self, point):
        # 时间模块提供零起始的月中心索引；文件名中的月份从 1 开始编号，故读取时加 1。
        # 跨年处仍用所选 MY 的首尾月循环，不读取或假装拥有相邻火星年的真实资料。
        season = point.seasonal_weights
        lower = self.load_month(season.lower_index + 1)
        upper = self.load_month(season.upper_index + 1)
        # 计算真太阳时在源 Time 数组里的左右采样及权重，24 时与 0 时按日周期衔接。
        diurnal = diurnal_weights(point.true_solar_time, self.hours)
        result = {}
        for name in lower:
            first, second = lower[name], upper[name]
            # 普通月平均量含日变化维：在两个相邻月份中各做一次线性日内插值。
            # RMS 源场没有这一日内维，保留其空间数组，仅参与后面的季节插值。
            if not name.startswith("rms_"):
                first = (first[diurnal.lower_index] * diurnal.lower_weight
                    + first[diurnal.upper_index] * diurnal.upper_weight)
                second = (second[diurnal.lower_index] * diurnal.lower_weight
                    + second[diurnal.upper_index] * diurnal.upper_weight)
            # 用月中心季节权重混合两月；结果去掉了 Time 维，但仍保留原水平和垂直索引。
            result[name] = first * season.lower_weight + second * season.upper_weight
        # 返回变量名到数组的字典；没有减去 RMS，也没有把插值背景标为每日真值。
        return result
