import math
import time
from datetime import datetime, timezone

# === Constants ===
class Constants:
    degrees_longitude_per_hour = 15
    degrees_to_radians = math.pi / 180
    earth_tilt_in_radians = 23.43715 * math.pi / 180
    julian_date_epoch = 2440587.5
    number_of_days_in_century = 36525
    number_of_seconds_per_day = 86400
    radians_to_degrees = 180 / math.pi
    three_sixty_degrees = 360


# === Coordinate Formatting ===
def decimal_to_deg_min_sec(lat, lon):
    def format_coord(value, is_lat):
        total_seconds = int(abs(value) * 3600)
        degrees = total_seconds // 3600
        minutes = (total_seconds % 3600) // 60
        direction = (
            "North" if (value >= 0 and is_lat) else
            "South" if (value < 0 and is_lat) else
            "East" if value >= 0 else "West"
        )
        return degrees, minutes, direction

    lat_deg, lat_min, lat_dir = format_coord(lat, is_lat=True)
    lon_deg, lon_min, lon_dir = format_coord(lon, is_lat=False)

    return f"{lat_deg:3d}° {lat_min:02d}' {lat_dir}  {lon_deg:3d}° {lon_min:02d}' {lon_dir}"


# === Astronomical Calculations ===
def julian_date(date):
    return Constants.julian_date_epoch + date.timestamp() / Constants.number_of_seconds_per_day

def julian_century(date):
    return (julian_date(date) - 2451545.0) / Constants.number_of_days_in_century

def geometric_mean_longitude_sun(t):
    return (280.46646 + t * 36000.76983 + t * t * 0.0003032) % Constants.three_sixty_degrees

def mean_anomaly(t):
    return 357.52911 + t * 35999.05029 - t * 0.0001537

def orbit_eccentricity(t):
    return 0.016708634 - t * (0.000042037 + 0.0000001267 * t)

def sun_equation_of_center(t):
    m = mean_anomaly(t)
    m_rad = m * Constants.degrees_to_radians
    return (
        math.sin(m_rad) * (1.914602 - t * (0.004817 + 0.000014 * t)) +
        math.sin(2 * m_rad) * (0.019993 - 0.000101 * t) +
        math.sin(3 * m_rad) * 0.000289
    )

def equation_of_time(date):
    t = julian_century(date)
    L0 = geometric_mean_longitude_sun(t)
    M = mean_anomaly(t)
    e = orbit_eccentricity(t)
    obliquity = 23.439 - 0.0000004 * t  # mean obliquity
    epsilon = obliquity * Constants.degrees_to_radians
    L0_rad = L0 * Constants.degrees_to_radians
    M_rad = M * Constants.degrees_to_radians

    y = math.tan(epsilon / 2)**2
    sin2L0 = math.sin(2 * L0_rad)
    sinM = math.sin(M_rad)
    cos2L0 = math.cos(2 * L0_rad)
    sin4L0 = math.sin(4 * L0_rad)
    sin2M = math.sin(2 * M_rad)

    E = y * sin2L0 - 2 * e * sinM + 4 * e * y * sinM * cos2L0 \
        - 0.5 * y * y * sin4L0 - 1.25 * e * e * sin2M

    return E * Constants.radians_to_degrees * 4  # in minutes


def subsolar_latitude(date):
    t = julian_century(date)
    L = geometric_mean_longitude_sun(t) + sun_equation_of_center(t)
    L_rad = L * Constants.degrees_to_radians
    lat_rad = math.asin(math.sin(L_rad) * math.sin(Constants.earth_tilt_in_radians))
    return lat_rad * Constants.radians_to_degrees

def subsolar_longitude(date):
    utc = date.astimezone(timezone.utc)
    utc_hours = utc.hour + utc.minute / 60 + utc.second / 3600
    eot = equation_of_time(date) / 60  # convert to hours
    lng = (12 - utc_hours - eot) * Constants.degrees_longitude_per_hour
    if lng > 180:
        lng -= 360
    elif lng < -180:
        lng += 360
    return lng


# === Main Loop ===
if __name__ == "__main__":
    delay_seconds = 5
    try:
        while True:
            now = datetime.now()
            lat = subsolar_latitude(now)
            lon = subsolar_longitude(now)
            coord_str = decimal_to_deg_min_sec(lat, lon)
            print(f"The subsolar point is now at: {coord_str} (at {now.strftime('%Y-%m-%d %H:%M:%S')})")
            time.sleep(delay_seconds)
    except KeyboardInterrupt:
        print("Loop interrupted by user.")
