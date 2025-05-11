import math
import time
from datetime import datetime, timezone, timedelta

class Constants:
    coordinates_string_format = "{:3d}° {:02d}' {}  {:3d}° {:02d}' {}"
    degrees_longitude_per_hour = 15
    degrees_to_radians = math.pi / 180
    earth_tilt_in_radians = 0.40905543478
    julian_date_for_jan_01_1970 = 2440587.5
    noon_time = 12
    number_of_days_in_a_century = 36525
    number_of_days_in_a_year = 365
    number_of_minutes_in_an_hour = 60
    number_of_minutes_in_a_day = 1440
    number_of_seconds_in_a_day = 86400
    number_of_seconds_in_a_year = number_of_seconds_in_a_day * number_of_days_in_a_year
    number_of_seconds_in_an_hour = 3600
    radians_to_degrees = 180 / math.pi
    three_sixty_degrees = 360
    one_eighty_degrees = 180
    number_of_hours_in_a_day = 24

def decimal_to_deg_min_sec(lat, lon):
    def convert(coord, pos_dir, neg_dir):
        total_seconds = int(coord * Constants.number_of_seconds_in_an_hour)
        degrees = total_seconds // Constants.number_of_seconds_in_an_hour
        seconds = abs(total_seconds % Constants.number_of_seconds_in_an_hour)
        minutes = seconds // Constants.number_of_minutes_in_an_hour
        direction = pos_dir if degrees >= 0 else neg_dir
        return abs(degrees), minutes, direction

    lat_d, lat_m, lat_dir = convert(lat, "North", "South")
    lon_d, lon_m, lon_dir = convert(lon, "East", "West")

    return Constants.coordinates_string_format.format(lat_d, lat_m, lat_dir, lon_d, lon_m, lon_dir)

def julian_date_from_datetime(dt):
    unix_time = dt.replace(tzinfo=timezone.utc).timestamp()
    return Constants.julian_date_for_jan_01_1970 + unix_time / Constants.number_of_seconds_in_a_day

def julian_century_since_2000(dt):
    return (julian_date_from_datetime(dt) - 2451545.0) / Constants.number_of_days_in_a_century

def geometric_mean_longitude_sun(t):
    return (280.46646 + t * 36000.76983 + t * t * 0.0003032) % Constants.three_sixty_degrees

def mean_anomaly(t):
    return 357.52911 + t * (35999.05029 - 0.0001537)

def eccentricity_of_earth_orbit(t):
    return 0.016708634 - t * (0.000042037 + 0.0000001267 * t)

def sun_equation_of_center(t):
    m = mean_anomaly(t)
    m_rad = m * Constants.degrees_to_radians
    return (math.sin(m_rad) * (1.914602 - t * (0.004817 + 0.000014 * t)) +
            math.sin(2 * m_rad) * (0.019993 - 0.000101 * t) +
            math.sin(3 * m_rad) * 0.000289)

def latitude_of_sun(dt):
    t = julian_century_since_2000(dt)
    geom_mean_long = geometric_mean_longitude_sun(t)
    true_long = geom_mean_long + sun_equation_of_center(t)
    true_long_rad = true_long * Constants.degrees_to_radians
    lat_rad = math.asin(math.sin(true_long_rad) * math.sin(Constants.earth_tilt_in_radians))
    return lat_rad * Constants.radians_to_degrees

def equation_of_time(dt):
    t = julian_century_since_2000(dt)
    L0 = geometric_mean_longitude_sun(t) * Constants.degrees_to_radians
    e = eccentricity_of_earth_orbit(t)
    M = mean_anomaly(t) * Constants.degrees_to_radians
    obliquity = 0.0430264916545165  # ~23.44 degrees in radians

    term1 = obliquity * math.sin(2 * L0)
    term2 = 2 * e * math.sin(M)
    term3 = 4 * e * obliquity * math.sin(M) * math.cos(2 * L0)
    term4 = 0.5 * obliquity * obliquity * math.sin(4 * L0)
    term5 = 1.25 * e * e * math.sin(2 * M)

    eot = 4 * (term1 - term2 + term3 - term4 - term5) * Constants.radians_to_degrees
    return eot

def subsolar_longitude(dt):
    eot = equation_of_time(dt)
    local_time = dt.astimezone()
    hour = local_time.hour + local_time.minute / 60.0
    seconds_from_gmt = local_time.utcoffset().total_seconds()

    if seconds_from_gmt <= 0:
        time_correction = 1
        day_correction = 0
    else:
        time_correction = -1
        day_correction = -Constants.number_of_hours_in_a_day

    gmt_hour = (hour - seconds_from_gmt / Constants.number_of_seconds_in_an_hour -
                time_correction * Constants.number_of_hours_in_a_day) % Constants.number_of_hours_in_a_day

    noon_delta = min(Constants.noon_time - gmt_hour - eot / Constants.number_of_minutes_in_an_hour + day_correction, 12.0)
    subsolar_lon = noon_delta * Constants.degrees_longitude_per_hour

    if subsolar_lon < -Constants.one_eighty_degrees and gmt_hour <= Constants.noon_time:
        lon_correction = Constants.one_eighty_degrees
    elif subsolar_lon <= -Constants.one_eighty_degrees and gmt_hour >= Constants.noon_time:
        lon_correction = 0
    elif gmt_hour >= Constants.number_of_hours_in_a_day:
        lon_correction = Constants.one_eighty_degrees
    else:
        lon_correction = 0

    return subsolar_lon + lon_correction

def get_subsolar_coordinates():
    now = datetime.now()
    lat = latitude_of_sun(now)
    lon = subsolar_longitude(now)
    return lat, lon

# Main loop
while True:
    now = datetime.now()
    lat, lon = get_subsolar_coordinates()
    coords = decimal_to_deg_min_sec(lat, lon)
    print(f"The subsolar point is now at: {coords} (at {now})")
    time.sleep(5)
