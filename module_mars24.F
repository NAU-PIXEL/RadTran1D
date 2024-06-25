module module_mars24

implicit none

real(kind(0.d0)) :: j2000_epoch = 2451545.0d0
real(kind(0.d0)) :: dtr = 0.017453292519943295d0

public ::&
     mars24_utc_to_tt_offset, mars24_julian_tt, mars24_j2000_offset_tt, &
     mars24_Mars_Mean_Anomaly, mars24_FMS_Angle, mars24_alpha_perturbs, &
     mars24_Equation_Of_Center, mars24_Mars_Ls, mars24_Equation_Of_Time, &
     mars24_j2000_ott_from_Mars_Solar_Date, mars24_Mars_Solar_Date, &
     mars24_Coordinated_Mars_Time, mars24_Local_Mean_Solar_Time, &
     mars24_Local_True_Solar_Time, mars24_subsolar_longitude, &
     mars24_solar_declination, mars24_heliocentric_distance, &
     mars24_heliocentric_longitude, mars24_heliocentric_latitude, &
     mars24_hourangle, mars24_solar_zenith, mars24_solar_elevation, &
     mars24_solar_azimuth, mars24_Mars_Year


contains

!====================================================================

  subroutine mars24_utc_to_tt_offset(jday, utc_to_tt_offset)
    !Returns the offset in seconds from a julian date in Terrestrial Time (TT)
    !to a Julian day in Coordinated Universal Time (UTC)
    real(kind(0.d0)) :: jday
    real(kind(0.d0)) :: off, utc_to_tt_offset
    off = 32.184d0 !!<1972
    
    if (jday .gt. 2441317.5d0) off=off+10d0 !1972 1 1
    if (jday .gt. 2441499.5d0) off=off+1.0d0 !1972 6 1
    if (jday .gt. 2441683.5d0) off=off+1.0d0 !1973 1 1
    if (jday .gt. 2442048.5d0) off=off+1.0d0 !1974 1 1
    if (jday .gt. 2442413.5d0) off=off+1.0d0 !1975 1 1
    if (jday .gt. 2442778.5d0) off=off+1.0d0 !1976 1 1
    if (jday .gt. 2443144.5d0) off=off+1.0d0 !1977 1 1
    if (jday .gt. 2443509.5d0) off=off+1.0d0 !1978 1 1
    if (jday .gt. 2443874.5d0) off=off+1.0d0 !1979 1 1
    if (jday .gt. 2444239.5d0) off=off+1.0d0 !1980 1 1
    if (jday .gt. 2444786.5d0) off=off+1.0d0 !1981 6 1
    if (jday .gt. 2445151.5d0) off=off+1.0d0 !1982 6 1
    if (jday .gt. 2445516.5d0) off=off+1.0d0 !1983 6 1
    if (jday .gt. 2446247.5d0) off=off+1.0d0 !1985 6 1
    if (jday .gt. 2447161.5d0) off=off+1.0d0 !1988 1 1
    if (jday .gt. 2447892.5d0) off=off+1.0d0 !1990 1 1
    if (jday .gt. 2448257.5d0) off=off+1.0d0 !1991 1 1
    if (jday .gt. 2448804.5d0) off=off+1.0d0 !1992 6 1
    if (jday .gt. 2449169.5d0) off=off+1.0d0 !1993 6 1
    if (jday .gt. 2449534.5d0) off=off+1.0d0 !1994 6 1
    if (jday .gt. 2450083.5d0) off=off+1.0d0 !1996 1 1
    if (jday .gt. 2450630.5d0) off=off+1.0d0 !1997 6 1
    if (jday .gt. 2451179.5d0) off=off+1.0d0 !1999 1 1
    if (jday .gt. 2453736.5d0) off=off+1.0d0 !2006 1 1
    if (jday .gt. 2454832.5d0) off=off+1.0d0 !2009 1 1
    if (jday .gt. 2456109.5d0) off=off+1.0d0 !2012 6 1

    !off= 66.184
    utc_to_tt_offset = off

  end subroutine mars24_utc_to_tt_offset

  subroutine mars24_julian_tt(jday_utc, julian_tt)
    !Returns the TT Julian day given a UTC Julian day
    real(kind(0.d0)) :: jday_utc
    real(kind(0.d0)) :: julian_tt, utc
    call mars24_utc_to_tt_offset(jday_utc, julian_tt)
    julian_tt = jday_utc + julian_tt/86400d0
  end subroutine mars24_julian_tt

  subroutine mars24_j2000_offset_tt(jday_tt, j2000_offset_tt)
    !Returns the julian day offset since the J2000 epoch
    real(kind(0.d0)) :: jday_tt
    real(kind(0.d0)) :: j2000_offset_tt
    j2000_offset_tt = jday_tt - j2000_epoch
  end subroutine mars24_j2000_offset_tt

  subroutine mars24_Mars_Mean_Anomaly(j2000_ott, Mars_Mean_Anomaly)
    !Calculates the Mars Mean Anomaly given a j2000 julian day offset
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: Mars_Mean_Anomaly
    !  real(kind(0.d0)) :: Mars_Mean_Anomaly

    Mars_Mean_Anomaly = 19.3870d0 + 0.52402075d0 * j2000_ott
    Mars_Mean_Anomaly = modulo(Mars_Mean_Anomaly+360d0,360d0)
  end subroutine mars24_Mars_Mean_Anomaly

  subroutine mars24_FMS_Angle(j2000_ott, FMS_Angle)
    !Returns the Fictional Mean Sun angle
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: FMS_Angle
    !  real(kind(0.d0)) :: FMS_Angle

    FMS_Angle = 270.3863d0 +  0.52403840d0 * j2000_ott
    FMS_Angle = modulo(FMS_Angle+360d0,360d0)
  end subroutine mars24_FMS_Angle

  subroutine mars24_alpha_perturbs(j2000_ott, alpha_perturbs)
    !Returns the perturbations to apply to the FMS Angle from orbital perturbations
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: alpha_perturbs

    real(kind(0.d0)), dimension(7) :: A, tau, phi

    integer :: i

    data A/0.0071d0, 0.0057d0, 0.0039d0, 0.0037d0, 0.0021d0, 0.0020d0, 0.0018d0/
    data tau/2.2353d0, 2.7543d0, 1.1177d0, 15.7866d0, 2.1354d0, 2.4694d0, 32.8493d0/
    data phi/49.409d0, 168.173d0, 191.837d0, 21.736d0, 15.704d0, 95.528d0, 49.095d0/

    alpha_perturbs =0.0d0

    do i=1,7 
       alpha_perturbs = alpha_perturbs + A(i) * cos(((0.985626d0*j2000_ott/tau(i)) + phi(i))*dtr)
    end do

  end subroutine mars24_alpha_perturbs

  subroutine mars24_Equation_Of_Center(j2000_ott, Equation_Of_Center)
    !The true anomaly (v) - the Mean anomaly (M)
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: M, pbs
    real(kind(0.d0)) :: Equation_Of_Center
    call mars24_Mars_Mean_Anomaly(j2000_ott, M)
    M=M*dtr
    call mars24_alpha_perturbs(j2000_ott, pbs)

    Equation_Of_Center = (10.691d0 + 3.0d-7 * j2000_ott)*sin(M) &
         + 0.6230d0 * sin(2*M) &
         + 0.0500d0 * sin(3*M) & 
         + 0.0050d0 * sin(4*M) & 
         + 0.0005d0 * sin(5*M) &
         + pbs

  end subroutine mars24_Equation_Of_Center

  subroutine mars24_Mars_Ls(j2000_ott, Mars_Ls)
    !Returns the Areocentric solar longitude
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: alpha, v_m
    real(kind(0.d0)) :: Mars_Ls
    
    call mars24_FMS_Angle(j2000_ott, alpha)
    call mars24_Equation_Of_Center(j2000_ott, v_m)

    Mars_Ls = (alpha+v_m)
    Mars_Ls = modulo(Mars_Ls+360d0, 360d0)

  end subroutine mars24_Mars_Ls

  subroutine mars24_Equation_Of_Time(j2000_ott, Equation_Of_Time)
    !Equation of Time
    real(kind(0.d0)) :: j2000_ott
    real(kind(0.d0)) :: Ls
    real(kind(0.d0)) :: Equation_Of_Time, eoc

    call mars24_Mars_ls(j2000_ott, Ls)
    Ls = Ls*dtr
    call mars24_Equation_Of_Center(j2000_ott, eoc)

    Equation_Of_Time = 2.861d0*sin(2*Ls) &
         - 0.071d0 * sin(4*ls) &
         + 0.002d0 * sin(6*ls) - eoc
  end subroutine mars24_Equation_Of_Time

  subroutine mars24_j2000_ott_from_Mars_Solar_Date(msd, j2_ott)
    real(kind(0.d0)) :: msd
    real(kind(0.d0)) :: j2, j2_ott
    j2 = ((msd*1.0d0 - 44796.0d0 + 0.00096d0)*1.027491252d0 + 4.5)
    call mars24_julian_tt(j2+j2000_epoch, j2_ott)
    j2_ott = j2_ott - j2000_epoch

  end subroutine mars24_j2000_ott_from_Mars_Solar_Date

  subroutine mars24_Mars_Solar_Date(j2000_ott, msd)
    real(kind(0.d0)) :: j2000_ott, msd
    msd = (((j2000_ott - 4.5)/1.027491252d0) + 44796.0d0 - 0.00096d0)
  end subroutine mars24_Mars_Solar_Date

  subroutine mars24_Coordinated_Mars_Time(j2000_ott, Coordinated_Mars_Time)
    !The Mean Solar Time at the Prime Meridian
    real(kind(0.d0)) :: j2000_ott, Coordinated_Mars_Time
    
    Coordinated_Mars_Time = 24.0d0 * (((j2000_ott - 4.5d0)/1.027491252d0) + 44796.0d0 - 0.00096d0)
    Coordinated_Mars_Time= modulo(24d0+Coordinated_Mars_Time, 24.0d0)

  end subroutine mars24_Coordinated_Mars_Time

  subroutine mars24_Local_Mean_Solar_Time(longitude, j2000_ott, Local_Mean_Solar_Time)
    !The Local Mean Solar Time given a planetographic longitude
    real(kind(0.d0)) :: longitude
    real(kind(0.d0)) :: j2000_ott, west_longitude
    real(kind(0.d0)) :: mtc, Local_Mean_Solar_Time
    call mars24_Coordinated_Mars_Time(j2000_ott, mtc)
    west_longitude = 360.0 - longitude !THIS OFFSET CHANGES TES LONGITUDES TO MOLA LONGITUDES +0.271d0
    Local_Mean_Solar_Time = mtc - longitude * (24d0/360d0)
    Local_Mean_Solar_Time = modulo(Local_Mean_Solar_Time + 24.0d0, 24.d0)

  end subroutine mars24_Local_Mean_Solar_Time


  subroutine mars24_Local_True_Solar_Time(longitude, j2000_ott, Local_True_Solar_Time)
    !Local true solar time is the Mean solar time + equation of time perturbation
    real(kind(0.d0)) :: j2000_ott, longitude
    real(kind(0.d0)) :: eot, lmst, Local_True_Solar_Time

    call mars24_Equation_Of_Time(j2000_ott, eot)
    call mars24_Local_Mean_Solar_Time(longitude, j2000_ott, lmst)

    Local_True_Solar_Time = lmst + eot*(24d0/360d0)
    Local_True_Solar_Time = modulo(Local_True_Solar_Time + 24d0, 24.0d0)

  end subroutine mars24_Local_True_Solar_Time

  subroutine mars24_subsolar_longitude(j2000_ott, subsolar_longitude)
    !    """returns the longitude of the subsolar point for a given julian day."""
    real(kind(0.d0)) :: j2000_ott, subsolar_longitude
    real(kind(0.d0)) :: MTC, EOT

    call mars24_Coordinated_Mars_Time(j2000_ott, mtc)
    call mars24_equation_of_time(j2000_ott, eot)
    eot = eot*24d0/360d0
    subsolar_longitude = (MTC + EOT)*(360d0/24d0) + 180d0
    subsolar_longitude = modulo(subsolar_longitude+360d0,360.d0)
  end subroutine mars24_subsolar_longitude

  subroutine mars24_solar_declination(ls, solar_declination)
    !    """Returns the solar declination"""
    real(kind(0.d0)) :: ls, solar_declination
    solar_declination = (asin(0.42565d0 * sin(ls*dtr)) + 0.25d0*dtr * sin(ls*dtr)) / dtr
  end subroutine mars24_solar_declination

  subroutine mars24_heliocentric_distance(j2000_ott, heliocentric_distance)
    !    """Instantaneous orbital radius"""
    real(kind(0.d0)) :: j2000_ott, heliocentric_distance
    real(kind(0.d0)) :: M
    call mars24_Mars_Mean_Anomaly(j2000_ott, M)
    M=M*dtr

    heliocentric_distance = 1.523679d0 * &
         (1.00436d0 - 0.09309d0*cos(M) &
         - 0.004336d0*cos(2*M) &
         - 0.00031d0*cos(3*M) &
         - 0.00003d0*cos(4*M))
  end subroutine mars24_heliocentric_distance

  subroutine mars24_heliocentric_longitude(j2000_ott, heliocentric_longitude)
    !    """Heliocentric longitude"""
    real(kind(0.d0)) :: j2000_ott, heliocentric_longitude
    real(kind(0.d0)) :: ls

    call mars24_Mars_Ls(j2000_ott, ls)

    heliocentric_longitude = ls + 85.061d0 - &
         0.015d0 * sin((71+2*ls)*dtr) - &
         5.5d-6*j2000_ott

    heliocentric_longitude = modulo(heliocentric_longitude + 360d0, 360d0)

  end subroutine mars24_heliocentric_longitude


  subroutine mars24_heliocentric_latitude(j2000_ott, heliocentric_latitude)
    !    """Heliocentric Latitude"""
    real(kind(0.d0)) :: j2000_ott, heliocentric_latitude
    real(kind(0.d0)) :: ls

    call mars24_Mars_Ls(j2000_ott, ls)


    heliocentric_latitude = -(1.8497d0 - 2.23d-5*j2000_ott) &
         * sin((ls - 144.50d0 + 2.57d-6*j2000_ott)*dtr)

  end subroutine mars24_heliocentric_latitude

  subroutine mars24_hourangle(longitude, j2000_ott, hourangle)
    !    """Hourangle is the longitude - subsolar longitude"""
    real(kind(0.d0)) :: j2000_ott, longitude, hourangle
    real(kind(0.d0)) :: subsolar_longitude

    call mars24_subsolar_longitude(j2000_ott, subsolar_longitude)
    hourangle = longitude*dtr - subsolar_longitude*dtr

  end subroutine mars24_hourangle

  subroutine mars24_solar_zenith(longitude,latitude, j2000_ott, solar_zenith)
    !    """Zenith Angle"""

    real(kind(0.d0)) :: j2000_ott, longitude, latitude
    real(kind(0.d0)) :: ha, ls, dec, cosZ, solar_zenith

    call mars24_hourangle(longitude, j2000_ott, ha)
    call mars24_Mars_Ls(j2000_ott, ls)
    call mars24_solar_declination(ls,dec)
    dec=dec*dtr

    cosZ = sin(dec) * sin(latitude*dtr) + &
         cos(dec)*cos(latitude*dtr)*cos(ha)

    solar_zenith = acos(cosZ)/dtr
  end subroutine mars24_solar_zenith

  subroutine mars24_solar_elevation(longitude, latitude, j2000_ott, solar_elevation)
    !    """Elevation = 90-Zenith"""
    real(kind(0.d0)) :: longitude, latitude, j2000_ott, solar_elevation, solar_zenith
    call mars24_solar_zenith(longitude, latitude, j2000_ott, solar_zenith)
    solar_elevation = 90 - solar_zenith

    
  end subroutine mars24_solar_elevation

  subroutine mars24_solar_azimuth(longitude, latitude, j2000_ott, solar_azimuth)
    !    """Azimuth Angle"""
    real(kind(0.d0)) :: longitude, latitude, j2000_ott
    real(kind(0.d0)) :: ha, ls, dec, denom, num, solar_azimuth

    call mars24_hourangle(longitude, j2000_ott, ha)
    call mars24_Mars_Ls(j2000_ott, ls)
    call mars24_solar_declination(ls, dec)
    dec=dec*dtr
    denom = (cos(latitude*dtr)*tan(dec) &
         - sin(latitude)*cos(ha))

    num = sin(ha) 
    solar_azimuth = modulo(360d0+atan2(num,denom)/dtr,360d0)

  end subroutine mars24_solar_azimuth

  subroutine mars24_Mars_Year(j2000_ott, MarsYear, day_of_year)

    ! """Returns the Mars Year date based on the reference date 1955 April 11, 10:56:31 mtc after 
    ! finding the j2k offsets of the zeroes of the Mars_Ls function. """
    ! Piqueux et al., Icarus, 251: 332-338, 2015

    integer, parameter :: num_years = 79
    real(kind(0.d0)), intent(in ) :: j2000_ott
    real(kind(0.d0)), intent(out) :: MarsYear, day_of_year
    real(kind(0.d0)) :: msd, msd_start_of_year
    real(kind(0.d0)), dimension(num_years) :: jday_vals, year_vals, year_length
    integer :: i

    data jday_vals / -16336.044076, -15649.093471, -14962.0892946, -14275.0960023,      &
     -13588.1458658, -12901.1772635, -12214.2082215, -11527.2637345, -10840.2842249,    & 
     -10153.2828749, -9466.3114025, -8779.3356111, -8092.3607738, -7405.4236452,        &
     -6718.4615347, -6031.4574604, -5344.4876509, -4657.5318339, -3970.5474528,         &
     -3283.5848372, -2596.6329362, -1909.6426682, -1222.6617049, -535.7040268,          &
     151.2736522, 838.2369682, 1525.1834712, 2212.1799182, 2899.1848518, 3586.1403058,  &
     4273.1024234, 4960.0765368, 5647.0207838, 6333.986502, 7020.9875066, 7707.9629132, &
     8394.9318782, 9081.9102062, 9768.8526533, 10455.8028354, 11142.8050514,            &
     11829.7873254, 12516.7417734, 13203.725449, 13890.6991502, 14577.6484912,          &
     15264.6324865, 15951.6217969, 16638.5798914, 17325.5517216, 18012.5209097,         &
     18699.4628887, 19386.4443201, 20073.4534421, 20760.4152811, 21447.3696661,         &
     22134.3466251, 22821.2966642, 23508.2529432, 24195.2539572, 24882.2400506,         &
     25569.2081296, 26256.1902459, 26943.1429481, 27630.0847446, 28317.0793316,         &
     29004.0710936, 29691.0238241, 30377.9991486, 31064.9784277, 31751.9249377,         &
     32438.896907, 33125.8902412, 33812.8520242, 34499.8183442, 35186.7944595,          &
     35873.740573, 36560.7112423, 37247.7247318 /

    data year_vals / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
            20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, &
            39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, &
            58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, &
            77, 78, 79 /

    data year_length / 686.95252, 686.950605, 687.0041764, 686.9932923, 686.9501365,    &
            686.9686023, 686.969042, 686.944487, 686.9795096, 687.00135, 686.9714724,   &
            686.9757914, 686.9748373, 686.9371286, 686.9621105, 687.0040743,            &
            686.9698095, 686.955817, 686.9843811, 686.9626156, 686.951901, 686.990268,  &
            686.9809633, 686.9576781, 686.977679, 686.963316, 686.946503, 686.996447,   &
            687.0049336, 686.955454, 686.9621176, 686.9741134, 686.944247, 686.9657182, &
            687.0010046, 686.9754066, 686.968965, 686.978328, 686.9424471, 686.9501821, &
            687.002216, 686.982274, 686.954448, 686.9836756, 686.9737012, 686.949341,   &
            686.9839953, 686.9893104, 686.9580945, 686.9718302, 686.9691881, 686.941979,&
            686.9814314, 687.009122, 686.961839, 686.954385, 686.976959, 686.9500391,   &
            686.956279, 687.001014, 686.9860934, 686.968079, 686.9821163, 686.9527022,  &
            686.9417965, 686.994587, 686.991762, 686.9527305, 686.9753245, 686.9792791, &
            686.94651, 686.9719693, 686.9933342, 686.961783, 686.96632, 686.9761153,    &
            686.9461135, 686.9706693, 687.0134895 /

    day_of_year = -999.d0
    if      (j2000_ott .lt. jday_vals(1)) then
                MarsYear = 1.d0+REAL(INT((j2000_ott-jday_vals(1))/year_length(1)))
    else if (j2000_ott .gt. jday_vals(num_years)) then
                MarsYear = 1.d0+REAL(INT((j2000_ott-jday_vals(num_years))/year_length(num_years)))
    else
          do i=2,num_years
            if ((jday_vals(i-1) .le. j2000_ott) .and. &
               (jday_vals(i)   .gt. j2000_ott)) then
                  MarsYear = year_vals(i-1)
                  call mars24_Mars_Solar_Date(jday_vals(i-1), msd_start_of_year)
                  call mars24_Mars_Solar_Date(j2000_ott, msd)
                  day_of_year = msd-msd_start_of_year
            endif
          enddo
    endif

  end subroutine mars24_Mars_Year

end module module_mars24
