MODULE inputs
  IMPLICIT NONE

  ! base case for weak 20
  REAL, PARAMETER ::  ws10 = 5.75,                           &
                      iTsurf = 25.0,                         &
                      hflux = -5.0e-7,                    &
                      ihb = 30.0,                           &
                      cd_fac = 0.7,                         &
                      c_alk = 1.5,                          &
                      c1  =  7.56903,              & ! Carbon Dioxide
                      c2  =  1.67006e03,             & !Bicarbonate
                      c3  =  3.14655e02,             & !carbonate
                      c4  =  2.96936e02,             & !BOH3
                      c5  =  1.18909e02,             & !BOH4
                      c6  =  6.30928e-03,           & !Hydrogen
                      c7  =  9.60492,             & !OH
                      ustokes = 5.75,                       &
                      Rgas = 0.0083143


END MODULE inputs
