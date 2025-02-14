import math
from datetime import datetime, timedelta, timezone
import time

class Constants:
    coordinates_string_format = "{:3d}° {:02d}' {}  {:3d}° {:02d}' {}"
    degrees_longitude_per_hour = 15
    degrees_to_radians = 0.01745329252
    earth_tilt_in_degrees = 23.43715
    earth_tilt_in_radians = 0.40905543478
    julian_date_for_jan011970_at0000gmt = 2440587.5
    noon_time = 12
    number_of_days_in_a_century = 36525
    number_of_days_in_a_year = 365
    number_of_hours_in_a_day = 24
    number_of_minutes_in_a_day = 1440
    number_of_minutes_in_a_year = number_of_days_in_a_year * number_of_minutes_in_a_day
    number_of_minutes_in_an_hour = 60
    number_of_seconds_in_a_day = 86400
    number_of_seconds_in_a_minute = 60
    number_of_seconds_in_a_year = number_of_seconds_in_a_day * number_of_days_in_a_year
    number_of_seconds_in_an_hour = 3600
    one_eighty_degrees = 180
    radians_to_degrees = 57.2957795131
    three_sixty_degrees = 360


class CoordinateConversions:
    @staticmethod
    def decimal_coordinates_to_deg_min_sec(latitude, longitude, format_str):
        lat_seconds = int(latitude * Constants.number_of_seconds_in_an_hour)
        lat_degrees = lat_seconds // int(Constants.number_of_seconds_in_an_hour)
        lat_seconds = abs(lat_seconds % int(Constants.number_of_seconds_in_an_hour))
        lat_minutes = lat_seconds // int(Constants.number_of_seconds_in_a_minute)
        lat_seconds %= int(Constants.number_of_seconds_in_a_minute)

        lon_seconds = int(longitude * Constants.number_of_seconds_in_an_hour)
        lon_degrees = lon_seconds // int(Constants.number_of_seconds_in_an_hour)
        lon_seconds = abs(lon_seconds % int(Constants.number_of_seconds_in_an_hour))
        lon_minutes = lon_seconds // int(Constants.number_of_seconds_in_a_minute)
        lon_seconds %= int(Constants.number_of_seconds_in_a_minute)

        lat_direction = "North" if lat_degrees >= 0 else "South"
        lon_direction = "East" if lon_degrees >= 0 else "West"

        return format_str.format(abs(lat_degrees), lat_minutes, lat_direction,
                                 abs(lon_degrees), lon_minutes, lon_direction)

    @staticmethod
    def deg_min_sec_coordinates_to_decimal(degrees, minutes, seconds, direction):
        sign = -1.0 if direction in ["S", "W"] else 1.0
        return (degrees + (minutes + seconds / Constants.number_of_seconds_in_a_minute) /
                Constants.number_of_seconds_in_a_minute) * sign


class Astrocalculations:
    @staticmethod
    def jd_from_date(date):
        jd = Constants.julian_date_for_jan011970_at0000gmt + date.timestamp() / Constants.number_of_seconds_in_a_day
        return jd

    @staticmethod
    def date_from_jd(jd):
        seconds_since_epoch = (jd - Constants.julian_date_for_jan011970_at0000gmt) * Constants.number_of_seconds_in_a_day
        date = datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=seconds_since_epoch)
        return date

    @staticmethod
    def julian_century_since_jan2000(date):
        jc = (Astrocalculations.jd_from_date(date) - 2451545) / Constants.number_of_days_in_a_century
        return jc

    @staticmethod
    def orbit_eccentricity_of_earth(t):
        ecc = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
        return ecc

    @staticmethod
    def mean_anomaly(t):
        m = 357.52911 + t * 35999.05029 - t * 0.0001537
        return m

    @staticmethod
    def equation_of_time(date):
        t = Astrocalculations.julian_century_since_jan2000(date)
        mean_longitude_sun_radians = Astrocalculations.geometric_mean_longitude_of_sun(t) * Constants.degrees_to_radians
        mean_anomaly_sun_radians = Astrocalculations.mean_anomaly(t) * Constants.degrees_to_radians
        eccentricity = Astrocalculations.orbit_eccentricity_of_earth(t)
        obliquity = 0.0430264916545165

        term1 = obliquity * math.sin(2 * mean_longitude_sun_radians)
        term2 = 2 * eccentricity * math.sin(mean_anomaly_sun_radians)
        term3 = 4 * eccentricity * obliquity * math.sin(mean_anomaly_sun_radians) * math.cos(2 * mean_longitude_sun_radians)
        term4 = 0.5 * obliquity * obliquity * math.sin(4 * mean_longitude_sun_radians)
        term5 = 1.25 * eccentricity * eccentricity * math.sin(2 * mean_anomaly_sun_radians)

        equation_of_time = 4 * (term1 - term2 + term3 - term4 - term5) * Constants.radians_to_degrees
        return equation_of_time

    @staticmethod
    def geometric_mean_longitude_of_sun(t):
        sun_geom_mean_long = 280.46646 + t * 36000.76983 + t * t * 0.0003032
        sun_geom_mean_long = sun_geom_mean_long % Constants.three_sixty_degrees
        return sun_geom_mean_long

    @staticmethod
    def latitude_of_sun(date):
        t = Astrocalculations.julian_century_since_jan2000(date)
        geom_mean_longitude = Astrocalculations.geometric_mean_longitude_of_sun(t)
        sun_true_longitude = geom_mean_longitude + Astrocalculations.sun_equation_of_center(t)
        sun_latitude = math.asin(math.sin(sun_true_longitude * Constants.degrees_to_radians) * math.sin(Constants.earth_tilt_in_radians))
        sun_latitude_degrees = sun_latitude * Constants.radians_to_degrees
        return sun_latitude_degrees

    @staticmethod
    def sun_equation_of_center(t):
        m = Astrocalculations.mean_anomaly(t)
        s_eoc = (math.sin(m * Constants.degrees_to_radians) * (1.914602 - t * (0.004817 + 0.000014 * t)) +
                 math.sin(2 * m * Constants.degrees_to_radians) * (0.019993 - 0.000101 * t) +
                 math.sin(3 * m * Constants.degrees_to_radians) * 0.000289)

        return s_eoc

    @staticmethod
    def sub_solar_longitude_of_sun_at_current_time(date):
        eot = Astrocalculations.equation_of_time(date)
        local_time = date.hour + date.minute / Constants.number_of_minutes_in_an_hour
        seconds_from_gmt = date.utcoffset().total_seconds() if date.utcoffset() else 0

        GMT = (local_time - seconds_from_gmt / Constants.number_of_seconds_in_an_hour) % Constants.number_of_hours_in_a_day
        noon_delta = (Constants.noon_time - GMT - eot / Constants.number_of_minutes_in_an_hour) % Constants.number_of_hours_in_a_day
        sub_solar_long = noon_delta * Constants.degrees_longitude_per_hour

        if sub_solar_long < -Constants.one_eighty_degrees:
            sub_solar_long += Constants.one_eighty_degrees

        return sub_solar_long

    @staticmethod
    def get_sub_solar_coordinates():
        now = datetime.now(timezone.utc)
        lat = Astrocalculations.latitude_of_sun(now)
        lon = Astrocalculations.sub_solar_longitude_of_sun_at_current_time(now)
        print(lat,lon)
        return lat, lon


def main():
    # secs_to_delay = 5
    # while True:
        subsolar_point = Astrocalculations.get_sub_solar_coordinates()



if __name__ == "__main__":
    main()