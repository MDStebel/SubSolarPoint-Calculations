//: ## Playground - Compute the Real-Time Subsolar Point Coordinates
//: ###### Based on the algorithm used in ISS Real-Time Tracker 3D.
//: ###### Created by Michael Stebel on 02/26/2021.
//: ###### Copyright © 2021-2025 Michael Stebel Consulting, LLC. All rights reserved.
import Foundation
import PlaygroundSupport

/// Variables and constants used by this playground
struct Constants {
    static let coordinatesStringFormat          = "%3d° %02d' %@  %3d° %02d' %@"
    static let degreesLongitudePerHour: Double  = 15
    static let degreesToRadians: Double         = 0.01745329252
    static let earthTiltInDegrees: Double       = 23.43715
    static let earthTiltInRadians: Double       = 0.40905543478
    static let julianDateForJan011970At0000GMT  = 2440587.5
    static let noonTime: Double                 = 12
    static let numberOfDaysInACentury: Double   = 36525
    static let numberOfDaysInAYear: Double      = 365
    static let numberOfHoursInADay: Double      = 24
    static let numberOfMinutesInADay: Double    = 1440
    static let numberOfMinutesInAYear: Double   = numberOfDaysInAYear * numberOfMinutesInADay
    static let numberOfMinutesInAnHour: Double  = 60
    static let numberOfSecondsInADay: Double    = 86400
    static let numberOfSecondsInAMinute: Double = 60
    static let numberOfSecondsInAYear: Double   = numberOfSecondsInADay * numberOfDaysInAYear
    static let numberOfSecondsInAnHour: Double  = 3600
    static let oneEightyDegrees: Double         = 180
    static let radiansToDegrees: Double         = 57.2957795131
    static let threeSixtyDegrees: Double        = 360
}

/// Format converters for coordinates
struct CoordinateConversions {
    /// Convert coordinates from decimal to degrees, minutes, seconds, and direction
    ///
    /// This is a format conversion only.
    /// - Parameters:
    ///   - latitude: Latitude as a Double
    ///   - longitude: Longitude as a Double
    ///   - format: The format to use in the conversion as a String
    /// - Returns: The coordinates string in deg min sec format
    static func decimalCoordinatesToDegMinSec(latitude: Double, longitude: Double, format: String) -> String {
        var latSeconds  = Int(latitude * Double(Constants.numberOfSecondsInAnHour))
        let latDegrees  = latSeconds / Int(Constants.numberOfSecondsInAnHour)
        latSeconds      = abs(latSeconds % Int(Constants.numberOfSecondsInAnHour))
        let latMinutes  = latSeconds / Int(Constants.numberOfSecondsInAMinute)
        latSeconds %= Int(Constants.numberOfSecondsInAMinute)

        var lonSeconds = Int(longitude * Double(Constants.numberOfSecondsInAnHour))
        let lonDegrees = lonSeconds / Int(Constants.numberOfSecondsInAnHour)
        lonSeconds     = abs(lonSeconds % Int(Constants.numberOfSecondsInAnHour))
        let lonMinutes = lonSeconds / Int(Constants.numberOfSecondsInAMinute)
        lonSeconds %= Int(Constants.numberOfSecondsInAMinute)

        return String(format: format, abs(latDegrees), latMinutes, { latDegrees >= 0 ? "North" : "South" }(), abs(lonDegrees), lonMinutes, { lonDegrees >= 0 ? "East" : "West" }())
    }

    /// Convert coordinates from degrees, minutes, seconds, and direction to decimal
    ///
    /// This is a format conversion only.
    /// - Parameters:
    ///   - degrees: Degrees as a Double
    ///   - minutes: Minutes as a Double
    ///   - seconds: Seconds as a Double
    ///   - direction: Direction as a String (either "N" or "S")
    /// - Returns: Decimal representation of coordinates as a Double
    static func degMinSecCoordinatesToDecimal(degrees: Double, minutes: Double, seconds: Double, direction: String) -> Double {
        let sign = (direction == "S" || direction == "W") ? -1.0 : 1.0

        return (degrees + (minutes + seconds / Double(Constants.numberOfSecondsInAMinute)) / Double(Constants.numberOfSecondsInAMinute)) * sign
    }
}

struct Astrocalculations {
    /// Convert a given Gregorian date to a Julian date
    /// - Parameter date: Gregorian date to convert as a Date
    /// - Returns: Julian date as a Double
    static func jDFromDate(date: Date) -> Double {
        let jD = Constants.julianDateForJan011970At0000GMT + date.timeIntervalSince1970 / Double(Constants.numberOfSecondsInADay)
        
        return jD
    }
    
    /// Convert a given Julian date to a Gregorian date
    /// - Parameter jd: Julian date as a Double
    /// - Returns: Standard date as a NSDate
    static func dateFromJd(jd: Double) -> NSDate {
        let gD = NSDate(timeIntervalSince1970: (jd - Constants.julianDateForJan011970At0000GMT) * Double(Constants.numberOfSecondsInADay))
        
        return gD
    }
    
    /// Calculate the Julian century since Jan-1-2000
    /// - Parameter date: Gregorian date as a Date
    /// - Returns: Julian century as a Double
    static func julianCenturySinceJan2000(date: Date) -> Double {
        let jc = (jDFromDate(date: date) - 2451545) / Double(Constants.numberOfDaysInACentury)
        
        return jc
    }
    
    /// Calculate the orbital eccentricity of the Earth
    ///
    /// The orbital eccentricity of an astronomical object is a dimensionless parameter that determines the amount by which its orbit around another body deviates from a perfect circle.
    /// - Parameter t: Julian century as a Double
    /// - Returns: Eccentricity as a Double
    static func orbitEccentricityOfEarth(t: Double) -> Double {
        let ecc = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
        
        return ecc
    }
    
    /// Calculate the Mean Anomaly of the Sun for a given date
    ///
    /// The mean anomaly is the angle between lines drawn from the Sun to the perihelion and to a point moving in the orbit at a uniform rate corresponding to the period of revolution of the planet.
    /// If the orbit of the planet were a perfect circle, then the planet as seen from the Sun would move along its orbit at a fixed speed.
    /// Then it would be simple to calculate its position (and also the position of the Sun as seen from the planet).
    /// The position that the planet would have relative to its perihelion if the orbit of the planet were a circle is called the mean anomaly.
    /// - Parameter t: Julian century as a Double
    /// - Returns: The mean anomaly as a Double
    static func meanAnomaly(t: Double) -> Double {
        let m = 357.52911 + t * 35999.05029 - t * 0.0001537
        
        return m
    }
    
    /// Calculate the Equation of Time for a given date
    ///
    /// The Equation of Time (EOT) is a formula used to convert between solar time and clock time, accounting for the Earth’s elliptical orbit around the Sun and its axial tilt.
    /// Essentially, the Earth does not move perfectly smoothly in a perfectly circular orbit, so the EOT adjusts for that.
    /// - Parameter date: A date as a Date type
    /// - Returns: Equation of time in minutes as a Double
    static func equationOfTime(for date: Date) -> Double {
        let t = julianCenturySinceJan2000(date: date)
        let meanLongitudeSunRadians = geometricMeanLongitudeOfSunAtCurrentTime(t: t) * Constants.degreesToRadians
        let meanAnomalySunRadians = meanAnomaly(t: t) * Constants.degreesToRadians
        let eccentricity = orbitEccentricityOfEarth(t: t)
        let obliquity = 0.0430264916545165
        
        let term1 = obliquity * sin(2 * meanLongitudeSunRadians)
        let term2 = 2 * eccentricity * sin(meanAnomalySunRadians)
        let term3 = 4 * eccentricity * obliquity * sin(meanAnomalySunRadians) * cos(2 * meanLongitudeSunRadians)
        let term4 = 0.5 * obliquity * obliquity * sin(4 * meanLongitudeSunRadians)
        let term5 = 1.25 * eccentricity * eccentricity * sin(2 * meanAnomalySunRadians)
        
        let equationOfTime = 4 * (term1 - term2 + term3 - term4 - term5) * Constants.radiansToDegrees
        
        return equationOfTime
    }
    
    /// Calculate the Sun Equation of Center
    ///
    /// The orbits of the planets are not perfect circles, but rather ellipses, so the speed of the planet in its orbit varies, and therefore the apparent speed of the Sun along the ecliptic also varies throughout the planet's year.
    /// The true anomaly is the angular distance of the planet from the perihelion of the planet, as seen from the Sun. For a circular orbit, the mean anomaly and the true anomaly are the same.
    /// The difference between the true anomaly and the mean anomaly is called the Equation of Center.
    /// - Parameter t: Julian century as a Double
    /// - Returns: The Sun equation of Center in radians
    static func sunEquationOfCenter(t: Double) -> Double {
        let m = meanAnomaly(t: t)
        let sEOC = sin(m * Double(Constants.degreesToRadians)) * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin(2 * m * Double(Constants.degreesToRadians)) * (0.019993 - 0.000101 * t) + sin(3 * m * Double(Constants.degreesToRadians)) * 0.000289
        
        return sEOC
    }
    
    /// Calculate the Geometric Mean Longitude of the Sun
    ///
    /// The mean longitude of the Sun, corrected for the aberration of light. Mean longitude is the ecliptic longitude at which an orbiting body could be found if its orbit were circular and free of perturbations.
    /// While nominally a simple longitude, in practice the mean longitude does not correspond to any one physical angle.
    /// - Parameter t: Julian century as a Double
    /// - Returns: The geometric mean longitude in degrees as a Double
    static func geometricMeanLongitudeOfSunAtCurrentTime(t: Double) -> Double {
        let sunGeometricMeanLongitude = Double(280.46646 + t * 36000.76983 + t * t * 0.0003032).truncatingRemainder(dividingBy: Constants.threeSixtyDegrees)
        
        return sunGeometricMeanLongitude
    }
    
    /// Calculate the exact current latitude of the Sun for a given date and time
    ///
    /// Based on the geometric mean longitude of the Sun.
    /// - Parameter date: A date as a Date type
    /// - Returns: Latitude in degrees as a Double
    static func latitudeOfSun(for date: Date) -> Double {
        let jC = julianCenturySinceJan2000(date: date)
        
        let geomMeanLongitude = geometricMeanLongitudeOfSunAtCurrentTime(t: jC)
        let sunTrueLongitude = geomMeanLongitude + Double(sunEquationOfCenter(t: jC))
        let latitudeOfSun = asin(sin(sunTrueLongitude * Constants.degreesToRadians) * sin(Constants.earthTiltInRadians))
        let sunTrueLatitudeInDegrees = latitudeOfSun * Constants.radiansToDegrees
        
        return sunTrueLatitudeInDegrees
    }
    
    /// Calculate the exact current longitude of the Sun for a given date and time
    ///
    /// This is an original algorithm that I created to calculate the subsolar point longitude, based on a given date and the equation of time.
    /// - Parameter date: A date as a Date type
    /// - Returns: The subsolar longitude as a Double
    static func subSolarLongitudeAtCurrentTimeUTC(for date: Date) -> Double {
        let utcCalendar = Calendar(identifier: .gregorian)
        let utcComponents = utcCalendar.dateComponents(in: TimeZone(secondsFromGMT: 0)!, from: date)

        guard let hour = utcComponents.hour,
              let minute = utcComponents.minute,
              let second = utcComponents.second else {
            return 0.0
        }

        // Convert UTC time to decimal hours
        let utcDecimalHours = Double(hour) + Double(minute) / Constants.numberOfMinutesInAnHour + Double(second) / Constants.numberOfSecondsInAnHour

        // Get the Equation of Time in minutes
        let eotMinutes = equationOfTime(for: date)
        let eotHours = eotMinutes / Constants.numberOfMinutesInAnHour

        // Subsolar longitude with Equation of Time correction
        var longitude = (12.0 - utcDecimalHours - eotHours) * Constants.degreesLongitudePerHour

        // Normalize to [-180°, 180°]
        if longitude > 180 {
            longitude -= 360
        } else if longitude < -180 {
            longitude += 360
        }

        return longitude
    }
    
    /// Get the subsolar coordinates at the current date and time
    ///
    /// The subsolar point is the position on Earth where the Sun is at the zenith.
    /// - Returns: The subsolar coordinates as a tuple of Doubles
    static func getSubSolarCoordinates() -> (latitude: Float, longitude: Float) {
        let now = Date()
        let lat = Float(latitudeOfSun(for: now))
        let lon = Float(subSolarLongitudeAtCurrentTimeUTC(for: now))
                        
        return (lat, lon)
    }
}

//: ### Run it in a loop
let secsToDelay: UInt32 = 5

while true {
    let now           = Date()
    let subsolarPoint = Astrocalculations.getSubSolarCoordinates()
    let lat           = Double(subsolarPoint.latitude)
    let lon           = Double(subsolarPoint.longitude)
    let coordinates   = CoordinateConversions.decimalCoordinatesToDegMinSec(latitude: lat, longitude: lon, format: Constants.coordinatesStringFormat)
    print("The subsolar point is now at: \(coordinates) (at \(now))")
    sleep(secsToDelay)
}
