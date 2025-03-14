MODULE inputs
  IMPLICIT NONE

  ! base case for weak 20 "hurricane" conditions
  REAL, PARAMETER ::  ws10 = 12.0,                           & !wind speed
                      iTsurf = 13,                         & !surface temperature (C)
                      hflux = 0.0e-7,                    &
                      ihb = 30.0,                           & !mixing depth (sort of)
                      cd_fac = 0.1,                         &
                      c_alk = 0,                          &
                      c1  =  0,              &  !  carbon dioxide (umol/kg)
                      c2  =  0,             &! bicarbonate (umol/kg)
                      c3  =  0,             & ! carbonate (umol/kg)
                      c4  =  0,             & ! boric acid (umol/kg)
                      c5= 0.3,                         & !phytoplankton (umolN/L*m^3/kg) L/m^3=0.001
                      c6=0.2,                          & !zooplankton (umolN/L*m^3/kg)
                      c7=1.6,                         & !nitrogen (umolN/l*m^3/kg)
                      ustokes = 12.0,                  & !stokes drift only used in speed2stress becomes u_10 (follow)
                      Rgas = 0.0083143 !R constant for gas


END MODULE inputs
