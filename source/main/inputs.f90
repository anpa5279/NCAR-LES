MODULE inputs
  IMPLICIT NONE

  ! base case for weak 20
  REAL, PARAMETER ::  ws10 = 5.75,                           &
                      iTsurf = 17.0,                         &
                      hflux = 0.0e-7,                    &
                      ihb = 33.0,                           &
                      cd_fac = 0.7,                         &
                      c_alk = 1.5,                          &
                      c1  =  7.56903,              & ! carbon dioxide
                      c2  =  1.67006e03,             & !bicraboante
                      c3  =  3.14655e02,             & !carboante
                      c4  =  2.96936e02,             & !boric acid
                      c5  =  1.18909e02,             & ! tetrahydroborate
                      c6  =  6.30928e-03,           & ! hydrogen ion
                      c7  =  9.60492,             & ! hydroxyl
                      !c8  =  2.5,             & !phytoplankton
                      !c9  =  1.5,           & ! zooplankton
                      !c10  =  4.0,             & ! nutrients
                      ustokes = 5.75,                       &
                      Rgas = 0.0083143


END MODULE inputs
