//
//  AstroCalculations.swift
//  ISS Real-Time Tracker 3D
//
//  Created by Michael Stebel on 10/26/20.
//  Copyright © 2020-2021 Michael Stebel Consulting, LLC. All rights reserved.
//

import SceneKit


/// A set of functions that provide astronomical calculations for solar coordinates.
struct AstroCalculations {
    
    /// Convert a Gregorian date to a Julian date
    /// - Parameter date: Gregorian date to convert
    /// - Returns: Julian date as a Double
    static func julianDate(date : Date) -> Double {
        
        let julianDateForJan011970 = 2440587.5
        
        return julianDateForJan011970 + date.timeIntervalSince1970 / 86400
        
    }
    
    
    /// Calculate the Julian century from a Gregorian date
    /// - Parameter date: Gregorian date
    /// - Returns: Julian century as a Double
    static func julianCentury(date: Date) -> Double {
        
        let jc = (julianDate(date: date) - 2451545) / Double(Globals.numberOfDaysInCentury)
        
        return jc
        
    }
    
    
    /// Calculate the Sun equation of Center
    /// - Parameter t: Julian century
    /// - Returns: The Sun equation of Center in radians
    static func sunEquationOfCenter(t: Double) -> Double {

        let m = ((357.52911 + t * (35999.05029 - t * 0.0001537)))           // The Sun geometric mean anomaly in degrees
        let cInRadians = sin(m * Double(Globals.degreesToRadians)) * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin(2 * m * Double(Globals.degreesToRadians)) * (0.019993 - 0.000101 * t) + sin(3 * m * Double(Globals.degreesToRadians)) * 0.000289

        return cInRadians
        
    }
        
    
    /// Calculate the exact geometric mean longitude of the Sun at the current time
    ///
    /// The mean longitude of the Sun, corrected for the aberration of light. The mean anomaly of the Sun (actually, of the Earth in its orbit around the Sun, but it is convenient to pretend the Sun orbits the Earth). Put in the range 0° to 360° by adding or subtracting multiples of 360° as needed.
    /// - Returns: The geometric mean longitude in degrees as a Float
    static func getGeometricMeanLongitudeOfSunAtCurrentTime() -> Float {
     
        let now = Date()
        let jC = julianCentury(date: now)

        let sunGeometricMeanLongitude = Float((280.46646 + jC * (36000.76983 + jC * 0.0003032))).truncatingRemainder(dividingBy: 360)
        
        return Float(sunGeometricMeanLongitude)
        
    }
    
    
    /// Calculate the current latitude of the Sun
    ///
    /// Based on the geometric mean longitude of the Sun
    /// - Returns: Latitude in degrees as a Float
    static func getLatitudeOfSunAtCurrentTime() -> Float {
        
        let geomMeanLongitude = getGeometricMeanLongitudeOfSunAtCurrentTime()
        let latitudeOfSun = asin(sin(geomMeanLongitude * Globals.degreesToRadians) * sin(Globals.earthTiltInRadians)) * Globals.radiansToDegrees

        return latitudeOfSun
        
    }
    
    
    /// Calculate the longitude of the Sun's current subsolar point
    ///
    /// The subsolar point is the location on the Earth where the Sun is directly overhead.
    /// This is based on an empirical algorithm, which I developed.
    /// - Returns: The current subsolar longitude as a Float
    static func getSubSolarLongitudeOfSunAtCurrentTime() -> Float {
        
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
 
}
