MODULE inputs
  IMPLICIT NONE

  ! base case for weak 20
  REAL, PARAMETER ::  ws10 = 5.75,                           &
                      iTsurf = 25.0,                         &
                      hflux = 0.0e-7,                    &
                      ihb = 33.0,                           &
                      cd_fac = 0.7,                         &
                      c_alk = 1.5,                          &
                      c1  =  7.56903,              &
                      c2  =  1.67006e03,             &
                      c3  =  3.14655e02,             &
                      c4  =  2.96936e02,             &
                      c5  =  1.18909e02,             &
                      c6  =  6.30928e-03,           &
                      c7  =  9.60492,             &
                      c8  =  2.5,             &
                      c9  =  1.5,           &
                      c10  =  4.0,             &
                      c11=  0.0,             &
                      ustokes = 5.75,                       &
                      Rgas = 0.0083143


END MODULE inputs
