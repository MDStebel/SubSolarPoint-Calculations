//: Various Astronomical Computations to Determine the Subsolar point coordinates
import Foundation


/// Static struct used as a namespace to hold global variables and constants
struct Globals {

    static let coordinatesStringFormat                      = "%3d° %02d' %@  %3d° %02d' %@"
    static let degreesLongitudePerHour: Float               = 15
    static let degreesToRadians: Float                      = .pi / 180
    static let earthTiltInDegrees: Float                    = 23.43715
    static let earthTiltInRadians: Float                    = earthTiltInDegrees * degreesToRadians
    static let numberOfDaysInAYear: Float                   = 365
    static let numberOfDaysInCentury: Double                = 36525
    
    static let numberOfSecondsInADay: Float                 = numberOfSecondsInAnHour * numberOfHoursInADay
    static let numberOfHoursInADay: Float                   = 24
    static let numberOfSecondsInAnHour: Float               = 3600
    static let numberOfMinutesInAnHour: Float               = 60
    static let numberOfSecondsInAMinute : Float             = 60
    static let noonTime: Float                              = 12
    static let radiansToDegrees: Float                      = 1 / degreesToRadians

}

/// Format converters for coordinates
struct CoordinateConversions {
    
    /// Convert coordinates from decimal to degrees, minutes, seconds, and direction
    ///
    /// This is a format conversion.
    /// - Parameters:
    ///   - latitude: Latitude as a Double.
    ///   - longitude: Longitude as a Double.
    ///   - format: String containing the format to use in the conversion.
    /// - Returns: The coordinates string in deg min sec format.
    static func decimalCoordinatesToDegMinSec(latitude: Double, longitude: Double, format: String) -> String {
        
        var latSeconds  = Int(latitude * Double(Globals.numberOfSecondsInAnHour))
        let latDegrees  = latSeconds / Int(Globals.numberOfSecondsInAnHour)
        latSeconds      = abs(latSeconds % Int(Globals.numberOfSecondsInAnHour))
        let latMinutes  = latSeconds / Int(Globals.numberOfSecondsInAMinute)
        latSeconds %= Int(Globals.numberOfSecondsInAMinute)
        
        var longSeconds = Int(longitude * Double(Globals.numberOfSecondsInAnHour))
        let longDegrees = longSeconds / Int(Globals.numberOfSecondsInAnHour)
        longSeconds     = abs(longSeconds % 3600)
        let longMinutes = longSeconds / Int(Globals.numberOfSecondsInAMinute)
        longSeconds %= Int(Globals.numberOfSecondsInAMinute)
        
        return String(format: format, abs(latDegrees), latMinutes, {return latDegrees >= 0 ? "North" : "South"}(), abs(longDegrees), longMinutes, {return longDegrees >= 0 ? "East" : "West"}())
        
    }
    
    
    /// Convert coordinates from degrees, minutes, seconds, and direction to decimal
    /// - Parameters:
    ///   - degrees: Degrees as a Double.
    ///   - minutes: Minutes as a Double.
    ///   - seconds: Seconds as a Double.
    ///   - direction: Direction as a String (either "N" or "S").
    /// - Returns: Decimal representation of coordinates as a Double
    static func degMinSecCoordinatesToDecimal(degrees: Double, minutes: Double, seconds: Double, direction: String) -> Double {
        
        let sign = (direction == "S" || direction == "W") ? -1.0 : 1.0
        
        return (degrees + (minutes + seconds / 60.0) / 60.0) * sign
        
    }
    
    
}


/// Convert a Gregorian date to a Julian date
/// - Parameter date: Gregorian date to convert
/// - Returns: Julian date as a double
func jdFromDate(date: Date) -> Double {
    
    let JD_JAN_1_1970_0000GMT = 2440587.5
    let jD = JD_JAN_1_1970_0000GMT + date.timeIntervalSince1970 / 86400.0

    return jD

}

//
//func dateFromJd(jd : Double) -> NSDate {
//    let JD_JAN_1_1970_0000GMT = 2440587.5
//    return  NSDate(timeIntervalSince1970: (jd - JD_JAN_1_1970_0000GMT) * 86400)
//}


public func julianCenturySinceJan2000(date: Date) -> Double {
    
    let now = Date()
    let jc = (jdFromDate(date: now) - 2451545) / Globals.numberOfDaysInCentury
    
    return jc
    
}


//public func sunEquationOfCenter(t: Double) -> Double {
//
//    let m = ((357.52911 + t * (35999.05029 - t * 0.0001537))) // Sun geometric mean anomaly
//    let cInRadians = sin(m * Double(Globals.degreesToRadians)) * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin(2 * m * Double(Globals.degreesToRadians)) * (0.019993 - 0.000101 * t) + sin(3 * m * Double(Globals.degreesToRadians)) * 0.000289
//
//    return cInRadians
//
//}


/// Calculate the mean longitude of the Sun at the current time
/// - Returns: mean longitude as a Float
public func getGeometricMeanLongitudeOfSunAtCurrentTime() -> Float {
 
    let now = Date()

    let julianCentury = julianCenturySinceJan2000(date: now)

    let sunLongitude = Float((280.46646 + julianCentury * (36000.76983 + julianCentury * 0.0003032))).truncatingRemainder(dividingBy: 360.0)
    
    return Float(sunLongitude)
    
}


/// Calculate the exact current latitude of the Sun
///
/// Based on the geometric mean longitude of the Sun
/// - Returns: Latitude in degrees as a Float
func getLatitudeOfSunAtCurrentTime() -> Float {
    
    let geomMeanLongitude = getGeometricMeanLongitudeOfSunAtCurrentTime()
    let latitudeOfSun = asin(sin(geomMeanLongitude * Globals.degreesToRadians) * sin(Globals.earthTiltInDegrees * Globals.degreesToRadians)) * Globals.radiansToDegrees
    
    return latitudeOfSun
    
}


/// Compute the Sun's current subsolar point longitude
///
/// The subsolar point is the location on Earth where the Sun is directly overhead.
/// - Returns: The subsolar longitude as a Float
func subSolarLongitudeOfSunAtCurrentTime() -> Float {
    
    var timeCorrection: Float
    var dayCorrection: Float
    var lonCorrection: Float
    
    // This determines current local and GMT time
    let localMins = Float(Calendar.current.component(.minute, from: Date()))
    let localHour = Float(Calendar.current.component(.hour, from: Date())) + localMins / Globals.numberOfMinutesInAnHour
    let secondsFromGMT = Float(TimeZone.current.secondsFromGMT())
    
    // Correct for time and day relative to GMT and the International Date Line
    if secondsFromGMT < 0 {
        timeCorrection = 1
        dayCorrection = 0
    } else {
        timeCorrection = -1
        dayCorrection = -Globals.numberOfHoursInADay
    }
    
    // Calculate GMT
    let GMT = localHour - secondsFromGMT / Globals.numberOfSecondsInAnHour - timeCorrection *  Globals.numberOfHoursInADay.truncatingRemainder(dividingBy: Globals.numberOfHoursInADay)
    
    // Now, calculate the displacement between current GMT and noontime in hours
    let noonHourDisplacement = Globals.noonTime - GMT + dayCorrection
    
    // The subsolar longitude is the displacement in hours times the number of degrees per hour (360/24=15)
    let subSolarLon = noonHourDisplacement * Globals.degreesLongitudePerHour
    
    // Now, determine if we've crossed the international date line. If so, we need to add 180 degrees.
    if subSolarLon < -180 && GMT <= Globals.noonTime {
        lonCorrection = 180
    } else if subSolarLon < -180 && GMT >= Globals.noonTime {
        lonCorrection = localHour >= GMT ? 0 : 180
    } else if GMT >= Globals.numberOfHoursInADay {
        lonCorrection = 180
    } else {
        lonCorrection = 0
    }

    return subSolarLon.truncatingRemainder(dividingBy: 180) + lonCorrection
    
}


// Run it
let coordinates = CoordinateConversions.decimalCoordinatesToDegMinSec(latitude: Double(getLatitudeOfSunAtCurrentTime()), longitude: Double(subSolarLongitudeOfSunAtCurrentTime()), format: Globals.coordinatesStringFormat)
print("The subsoloar point is at: \(coordinates)")
