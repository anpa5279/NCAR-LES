MODULE inputs
  IMPLICIT NONE

  ! base case for weak 20 "hurricane" conditions
  REAL, PARAMETER ::  ws10 = 12.0,                           & !wind speed
                      iTsurf = 13,                         & !surface temperature (C)
                      hflux = 0.0e-7,                    &
                      ihb = 30.0,                           & !mixing depth (sort of)
                      cd_fac = 0.1,                         &
                      c_alk = 1.5,                          &
                      c1  =  7.56903,              &  ! here and below are all the organic compounds simulated
                      c2  =  1.67006e03,             &
                      c3  =  3.14655e02,             &
                      c4  =  2.96936e02,             &
                      c5  =  1.18909e02,             &
                      c6  =  6.30928e-03,           &
                      c7  =  9.60492,             &
                      ustokes = 12.0,                       & !stokes drift only used in speed2stress becomes u_10 (follow)
                      Rgas = 0.0083143 !R constant for gas


END MODULE inputs
