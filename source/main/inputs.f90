MODULE inputs
    IMPLICIT NONE

    ! base case for weak 20
    REAL, PARAMETER :: ws10 = 5.75, &
                       iTsurf = 25.0, &
                       hflux = 0.0e-7, &
                       ihb = 30.0, &
                       cd_fac = 0.7, &
                       c_alk = 1.5, &
                       c1 = 7.57, & ! carbon dioxide
                       c2 = 1.67e03, & ! bicarboante
                       c3 = 3.15e02, & ! carbonate
                       c4 = 2.97e02, & ! B(OH)3
                       c5 = 1.19e02, & ! B(OH)4
                       c6 = 6.30928e-03, & !
                       c7 = 9.6, & ! hydroxide 
                       ustokes = 5.75, &
                       Rgas = 0.0083143

END MODULE inputs
