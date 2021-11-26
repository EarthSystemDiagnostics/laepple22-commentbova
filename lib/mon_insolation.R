


##' @title Calculate 'monthly' ( = 12 roughly equally sized annual bits) and annual insolation
##' @param kyear vector of time (kyr BP)
##' @param lat scalar latitude in degN
##' @return vector with 13 values; 'month' 1:12 and annual in [13]
##' @author Thomas Laepple 
mon_insolation<-function(kyear,lat)
{
    ##Hardcoded days of each month; Bova uses a sliding 30day window; Details of the seasonal definition
    #are not relevant as also shown by Bova
    MONTHBREAKS<-c(0,31,59,90,120,151,181,212,243,273,304,334,365)
    
    if (sum(kyear<0)>0) {stop("negative age") 
        }
    
    daily<-matrix(NA,length(kyear),365)
    monthly<-matrix(NA,length(kyear),13)
    for (i in 1:length(kyear)) daily[i,]<-DailyInsolation(kyear[i],lat,1:365)$Fsw
   
    for (iMon in 1:12)
    {
        rangeDay<-(MONTHBREAKS[iMon]+1):(MONTHBREAKS[iMon+1])
        monthly[,iMon]<-rowMeans(daily[,rangeDay])
    }
    monthly[,13]<-rowMeans(daily)
    
    return(monthly)
}


#' @title Calculate the Eccentricity anomaly
#' @param M mean anomaly
#' @param e
#' @param tol tolerance
#' @description Calculate the Eccentricity anomaly E from the mean anomaly M, iteratively
#' inversion of Keplers formula: E-e*sin(E)=M
#' based on Grossman, N.: Orbits under the Inverse Square Law, in
#' The Sheer Joy of Celestial Mechanics, edited by N. Grossman, pp. 41–61,
#' Birkhäuser Boston, Boston, MA., 1996.
#' @return Eccentricity anomaly E
#' @author Thomas Laepple
calcE <- function(M, e, tol = 1e-8)
{
  E <- vector()
  for (i in 1:length(M))
  {
    Etemp <- M[i]
    ratio <- 1
    while (abs(ratio) > tol)
    {
      f.E <- Etemp - e * sin(Etemp) - M[i]
      f.Eprime <- 1 - e * cos(Etemp)
      ratio <- f.E / f.Eprime
      if (abs(ratio) > tol)
        Etemp <- Etemp - ratio
      else
        E[i] <- Etemp
    }
  }
  return(E)
}


#' @title Calculate the true anomaly
#' @param E
#' @param e
#' @description Calculate the true anomaly f (theta) from the eccentricity
#' anomaly true anomaly = angle between the point and perihelion
#' based on Grossman, N.: Orbits under the Inverse Square Law, in
#' The Sheer Joy of Celestial Mechanics, edited by N. Grossman, pp. 41–61,
#' Birkhäuser Boston, Boston, MA., 1996.
#' @return true anomaly f (theta)
#' @author Thomas Laepple
calcf <- function(E, e){
  temp <- atan(sqrt((1 + e) / (1 - e)) * tan(0.5 * E)) * 2
  temp[temp < 0] <- temp[temp < 0] + 2 * pi
  return(temp)
}


#' @title days to true angle
#' @param ecc eccentricity
#' @param f true angle from perihelion
#' @references Grossman, N.: Orbits under the Inverse Square Law, in
#' The Sheer Joy of Celestial Mechanics, edited by N. Grossman, pp. 41–61,
#' Birkhäuser Boston, Boston, MA., 1996.
#' @return days from perihelion, calculated using the mean angle M
#' @author Thomas Laepple
t.from.f <- function(ecc, f) {

  #length of tropical year in days
  Tyear <- 365.2422

  #Calculate the eccentricity anomaly E
  E <- atan(tan(0.5 * f) / sqrt((1 + ecc) / (1 - ecc))) * 2 #(181)

  #Keplers equation to get the mean anomaly (188)
  M <-  E - ecc * sin(E)

  #Convert into days
  n = 2 * pi / Tyear
  t <- M / n #182
  if (t < 0)
    t <- t + Tyear
  return(t)
}


#' @title true angle from days
#' @param ecc eccentricity
#' @param t days from perihelion
#' @return  true angle from perihelion (via the mean angle M)
#' @references Grossman, N.: Orbits under the Inverse Square Law, in
#' The Sheer Joy of Celestial Mechanics, edited by N. Grossman, pp. 41–61,
#' Birkhäuser Boston, Boston, MA., 1996.
#' @author Thomas Laepple
f.from.t <- function(ecc, t)
{
  #length of tropical year in days
  Tyear <- 365.2422
  n = 2 * pi / Tyear
  M <- t * n
  E <- calcE(M, ecc)
  f <- calcf(E, ecc)
  return(f)
}






#'@title  Computes daily average insolation as a function of day and latitude
#'  given the orbital parameters
#'@description Insolation calculation, converted and adapted from the Eisenmann and
#'  Huybers,2006 Matlab Code, based on Berger 1991
#'
#'@section day type details:
#'
#'day_type=1 (default): day input is calendar day (1-365.24), where day 1 is
#' January first. The calendar is referenced to the vernal equinox which always
#' occurs at day 80.
#'
#' day_type=2: day input is solar longitude (0-360 degrees). Solar longitude is
#' the angle of the Earth's orbit measured from spring equinox (21 March). Note
#' that calendar days and solar longitude are not linearly related because, by
#' Kepler's Second Law, Earth's angular velocity varies according to its
#' distance from the sun.
#'
#' day_type=3: Custom calendar that alpha (angle since equinox) is reached at
#' day T.alpha this uses an iterative solution of Keplers equation and is
#' therefore much slower than the day_type=1 which uses a trigeometric
#' approximation of Kepler's equation.
#'
#' @param lat Latitude in degrees (-90 to 90).
#' @param day  Indicator of time of year; calendar day by default.
#' @param ecc  eccentricity
#' @param obliquity obliquity
#' @param long_perh longitude of perihelion
#' @param day_type  Convention for specifying time of year (+/- 1,2) [optional], see details
#' @param alpha parameter for custom calendar,see details
#' @param T.alpha parameter for custom calendar,see details
#' @return list(Fsw=Fsw, ecc=ecc, obliquity=obliquity, long_perh=long_perh, lambda=lambda/2/pi*360)
#' Fsw = Daily average solar radiation in W/m^2.
#'   ... and all orbital parameters
#' @author Thomas Laepple
#' @export
DailyInsolationParam<-function(lat, day, ecc, obliquity, long_perh, day_type=1,
                               alpha=0, T.alpha=80)
{
  epsilon <- obliquity * pi/180
  omega <- long_perh * pi/180

  # Calculate insolation -----
  lat = lat * pi / 180 # latitude

  # lambda (or solar longitude) is the angular distance along Earth's orbit
  # measured from spring equinox (21 March)

  if (day_type == 1)
    # calendar days
  {
    # estimate lambda from calendar day using an approximation from Berger 1978 section 3
    delta_lambda_m <- (day - 80) * 2 * pi / 365.2422
    ## das waere lambda bei gleich langen Tagen

    beta <- (1 - ecc^2)^(1/2)

    lambda_m0 <- (-2) * ((1/2 * ecc + 1/8 * ecc^3) * (1 + beta) * sin(-omega) -
                           1/4 * ecc^2 * (1/2 + beta) * sin(-2 * omega) + 1/8 *
                           ecc^3 * (1/3 + beta) * (sin(-3 * omega)))

    lambda_m <- lambda_m0 + delta_lambda_m

    lambda <- lambda_m + (2 * ecc - 1/4 * ecc^3) * sin(lambda_m - omega) +
      (5/4) * ecc^2 * sin(2 * (lambda_m - omega)) + (13/12) * ecc^3 *
      sin(3 * (lambda_m - omega))


  } else if (day_type==2)
    #solar longitude (1-360)
  {lambda = day*2*pi/360
  # lambda=0 for spring equinox
  } else if (day_type==3) {
    #t.from.f(eccc,alpha-omega) = days since Perihelion in which
    #the angle alpha, measured from the equinox is reached
    #tau
    tau = -1*t.from.f(ecc,alpha-omega) + T.alpha
    #angle since equinox given the days with the calendar fixed using
    #the condition that alpha (angle since equinox) is reached at day T.alpha
    lambda<-f.from.t(ecc,day-tau)+omega} else {
      stop('Error: invalid day_type')
    }

  So <- 1365  # solar constant (W/m^2)
  delta <- asin(sin(epsilon) * sin(lambda))  # declination of the sun
  Ho <- acos(-tan(lat) * tan(delta))  # hour angle at sunrise/sunset

  # no sunrise or no sunset: Berger 1978 eqn (8),(9)

  Ho[(abs(lat) >= (pi/2 - abs(delta))) & (lat * delta > 0)] <- pi
  Ho[(abs(lat) >= (pi/2 - abs(delta))) & (lat * delta <= 0)] <- 0

  # Insolation: Berger 1978 eq (10)
  Fsw = So/pi * (1 + ecc * cos(lambda - omega))^2/(1 - ecc^2)^2 *
    (Ho * sin(lat) *  sin(delta) + cos(lat) * cos(delta) * sin(Ho))

  return(list(Fsw = Fsw, ecc = ecc, obliquity = obliquity,
              long_perh = long_perh, lambda = lambda/2/pi * 360))
}


#' @title Annual insolation for any year of the last 5Mio years
#' @param kyear Thousands of years before present (0 to 5000).
#' @param lat Latitude in degrees (-90 to 90).
#' @return daily average solar radiation in W/m^2 averaged over the year
#' @author Thomas Laepple
#' @details for details see daily_insolation_param
#' @examples
#' plot(0:1000, AnnualInsolation(0:1000, 40), main = "mean annual insolation at 40degN",
#' xlab = "kyr BP", ylab = "W/m2", type = "l")
#' @export
AnnualInsolation <- function(kyear, lat) {
  result <- vector()
  for (i in 1:length(kyear)) result[i] <- mean(DailyInsolation(kyear[i],
                                                               lat, 1:365)$Fsw)
  return(result)
}



#' @title Daily average insolation for any time in the last 5 Million years
#' @description Computes the daily average insolation, converted and adapted from
#' the Eisenmann and Huybers,2006 Matlab Code.
#'
#' This script contains orbital parameter data for the past 50000 years
#'   from Berger and Loutre (1991).
#'
#' Detailed description of calculation:
#'   Values for eccentricity, obliquity, and longitude of perihelion for the
#'   past 5 Myr are taken from Berger and Loutre 1991 (data from ncdc.noaa.gov).
#'   If using calendar days, solar longitude is found using an approximate
#'   solution to the differential equation representing conservation of angular
#'   momentum (Kepler's Second Law) (day_type=1), or an iterative solution of
#'   Keplers Equation (day_type=3).
#'
#'   Given the orbital parameters and solar longitude, daily average insolation
#'   is calculated exactly following Berger 1978.
#'
#' Orbital parameters from Berger and Loutre 1991.
#'
#' Calculation of daily average insolation following Berger A. (1978).
#'
#' Authors of original MATLAB version :
#'   Ian Eisenman and Peter Huybers, Harvard University, August 2006
#'   eisenman@fas.harvard.edu
#'   This file is available online at
#'   http://deas.harvard.edu/~eisenman/downloads
#'
#' Translation into R and extension to arbitray calendars; Thomas Laepple, 2008
#'
#' @section day type details:
#'
#' day_type=1 (default): day input is calendar day (1-365.24), where day 1 is
#' January first. The calendar is referenced to the vernal equinox which always
#' occurs at day 80.
#'
#' day_type=2: day input is solar longitude (0-360 degrees). Solar longitude is
#' the angle of the Earth's orbit measured from spring equinox (21 March). Note
#' that calendar days and solar longitude are not linearly related because, by
#' Kepler's Second Law, Earth's angular velocity varies according to its
#' distance from the sun.
#'
#' day_type=3: Custom calendar day input is calendar day (1-365.24), where day 1
#' is January first. The calendar is referenced to the solar longitude alpha
#' (0..2pi, angle of the earth orbit from spring equinox) which always occurs at
#' day T.alpha. This uses an iterative solution of Keplers equation and is therefore
#' much slower than the day_type=1 which uses a trigeometric approximation
#' of Kepler's equation.
#'
#' @references Berger A. and Loutre M.F. (1991). Insolation values for the
#' climate of the last 10 million years. Quaternary Science Reviews, 10(4),
#' 297-317.
#'
#' Berger A. (1978). Long-term variations of daily insolation and Quaternary
#' climatic changes. Journal of Atmospheric Science, 35(12), 2362-2367.
#'
#' Grossman, N.: Orbits under the Inverse Square Law, in
#' The Sheer Joy of Celestial Mechanics, edited by N. Grossman, pp. 41–61,
#' Birkhäuser Boston, Boston, MA., 1996.
#'
#' @param kyear  Thousands of years before present (0 to 5000).
#' @param lat  Latitude in degrees (-90 to 90).
#' @param day  Indicator of time of year; calendar day by default.
#' @param day_type Convention for specifying time of year (1,2,3) [optional], see details
#' @param fast TRUE = use precalculated values of orbital parameters in 100yr steps
#' @param T.alpha for custom calendar (day_type=3)
#' @param alpha for custom calendar (day_type=3)
#' @return list(Fsw=Fsw,ecc=ecc,obliquity=obliquity,long_perh=long_perh,lambda=lambda/2/pi*360)
#'   Fsw = Daily average solar radiation in W/m^2.
#'   ... and all orbital parameters
#' @author Thomas Laepple
#' @examples
#' par(mfcol = c(2, 1))
#' plot(0:1000, DailyInsolation(0:1000, 40, 180)$Fsw, main = "day 180 (= peak summer) insolation at 40degN",
#'      xlab = "kyr BP", ylab = "W/m2", type = "l")
#' plot(1:365, DailyInsolation(0, 40, 1:365)$Fsw, main = "modern daily insolation at 40degN",
#'      xlab = "day of the year", ylab = "W/m2", type = "l")
#' @export
DailyInsolation<-function(kyear,lat,day,day_type=1,fast=TRUE,T.alpha=80,alpha=0) {

  # === Get orbital parameters ===

  if (fast) temp<-orbital_parameters_fast(kyear) else  temp<-orbital_parameters(kyear)
  ecc<-temp$ecc
  epsilon<-temp$epsilon
  omega<-temp$omega

  # For output of orbital parameters
  obliquity<-epsilon*180/pi
  long_perh<-omega*180/pi

  # === Calculate insolation ===
  lat=lat*pi/180 # latitude

  # lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
  if (day_type ==1)
    # calendar days
  {
    # estimate lambda from calendar day using an approximation from Berger 1978 section 3
    delta_lambda_m <- (day - 80) * 2 * pi/365.2422
    beta <- (1 - ecc^2)^(1/2)
    lambda_m0 <- (-2) * ((1/2 * ecc + 1/8 * ecc^3) * (1 + beta) * sin(-omega) - 1/4 *
                           ecc^2 * (1/2 + beta) * sin(-2 * omega) + 1/8 * ecc^3 * (1/3 + beta) * (sin(-3 *
                                                                                                        omega)))
    lambda_m <- lambda_m0 + delta_lambda_m
    lambda <- lambda_m + (2 * ecc - 1/4 * ecc^3) * sin(lambda_m - omega) + (5/4) * ecc^2 *
      sin(2 * (lambda_m - omega)) + (13/12) * ecc^3 * sin(3 * (lambda_m - omega))
  } else if (day_type==2)
    #solar longitude (1-360)
  {
    # lambda=0 for spring equinox
    lambda=day*2*pi/360
  } else if (day_type==3)
  {
    #t.from.f(eccc,alpha-omega) = days since perihelion in which
    #the angle alpha, measured from the equinox is reached

    #tau
    tau = -1*t.from.f(ecc,alpha-omega) + T.alpha

    #angle since equinox given the days with the calendar fixed using
    #the condition that alpha (angle since equinox) is reached at day T.alpha
    lambda<-f.from.t(ecc,day-tau)+omega



  }else stop('Error: invalid day_type')

  So<-1365 # solar constant (W/m^2)
  delta<-asin(sin(epsilon)*sin(lambda)) # declination of the sun
  Ho<-acos(-tan(lat)*tan(delta)) # hour angle at sunrise/sunset

  # no sunrise or no sunset: Berger 1978 eqn (8),(9)

  Ho[(abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta > 0 )]  <- pi
  Ho[( abs(lat) >= (pi/2 - abs(delta)) ) & ( lat*delta <= 0 )] <- 0

  # Insolation: Berger 1978 eq (10)
  Fsw=So/pi*(1+ecc*cos(lambda-omega))^2 /(1-ecc^2)^2 * ( Ho*sin(lat)*sin(delta) + cos(lat)*cos(delta)*sin(Ho))

  return(list(Fsw=Fsw,ecc=ecc,obliquity=obliquity,long_perh=long_perh,lambda=lambda/2/pi*360 ))
}




#' @title return orbital parameters for a given year using preloaded 100yr resolution values
#' @param kyear  time in 1000 years before present (0..5000)
#' @return list(ecc=ecc,epsilon=epsilon,omega=omega)
#' @author Thomas Laepple
orbital_parameters_fast<-function(kyear){
  if ( length(orbital_global$ecc) != 50001) stop("orbital parameters not initalized")

  return(
    list(ecc=orbital_global$ecc[kyear*10+1],
              epsilon=orbital_global$epsilon[kyear*10+1],
              omega=orbital_global$omega[kyear*10+1])
         )

}



#' @title Load the orbital parameters and return them for a given year
#' @details
#' Load orbital parameters (given each kyr for 0-5Mya)
#' Load the matrix contains data from Berger and Loutre (1991),
#'   downloaded as ORBIT91 from ncdc.noaa.gov
#' @param kyear time in 1000 years before present (0..5000)
#' @param FILEDATA Filename
#' @return list(ecc=ecc,epsilon=epsilon,omega=omega)
#' @author Thomas Laepple
orbital_parameters <- function(kyear, FILEDATA = "e:/data/insolation/data/ins_data.txt") {
  # kyears before present for data (kyear0>=0);
  kyear0 <- -1 * orbit91[, 1]
  # eccentricity
  ecc0 <- orbit91[, 2]

  # add 180 degrees to omega (see lambda definition, Berger 1978 Appendix)
  # longitude of perihelion (precession angle)
  omega0 <- orbit91[, 3] + 180
  # remove discontinuities (360 degree jumps)
  omega0 = unwrap(omega0 * pi/180) * 180/pi

  epsilon0 = orbit91[, 4]  # obliquity angle

  # Interpolate to requested dates
  ecc = spline(x = kyear0, y = ecc0, xout = kyear)$y
  omega = spline(x = kyear0, y = omega0, xout = kyear)$y * pi/180
  epsilon = spline(x = kyear0, y = epsilon0, xout = kyear)$y * pi/180

  return(list(ecc = ecc, epsilon = epsilon, omega = omega))
}


#' @title  Correct the radian phase angles
#' @details Q = unwrap(P) corrects the radian phase angles in array P by adding
#' multiples of \eqn{+- 2\pi} when absolute jumps between consecutive array
#' elements are greater than pi radians.
#' based on http://ccrma.stanford.edu/~jos/sasp/Matlab_listing_unwrap_m.html
#' @param p radian phase angles array
#' @return unwrapped angles
#' @author Thomas Laepple
unwrap <- function(p) {

  N = length(p)
  up = rep(0, N)
  pm1 = p[1]
  up[1] = pm1
  po = 0
  thr = pi
  pi2 = 2 * pi
  for (i in 2:N) {
    cp <- p[i] + po
    dp <- cp - pm1
    pm1 <- cp
    if (dp >= thr) {
      while (dp >= thr) {
        po <- po - pi2
        dp <- dp - pi2
      }
    }
    if (dp <= ((-1) * thr)) {
      while (dp <= thr) {
        po <- po + pi2
        dp <- dp + pi2
      }
    }
    cp = p[i] + po
    pm1 = cp
    up[i] = cp
  }
  return(up)
}

#Load the matrix contains data from Berger and Loutre (1991),
##'   downloaded as ORBIT91 from ncdc.noaa.gov
#orbit91<-read.table(paste(path,"ins_data.txt",sep=""))
#use_data(orbit91)
#usethis::use_data(orbital_global, internal = TRUE)
#orbital_global<-orbital_parameters((0:50000)/10)








