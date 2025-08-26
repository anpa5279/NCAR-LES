SUBROUTINE rhs_scl(istep, iscl)
! GET RHS OF SCALAR EQUATION (ISCL) MONOTONE SCALAR FLUXES ONLY IN Z FOR PENCIL
! SIZE (NNX,IYS:IYE,IZS:IZE)
! CARE IS TAKEN SO THAT IF MONOTONE IS ON THEN CONSERVATIVE HORIZONTAL FORM IS
! USED

    USE pars
    USE inputs
    USE fields
    USE fftwk
    USE con_data
    USE con_stats

    REAL :: fnt1(nnx, iys:iye, izs:ize)
    REAL :: tx(nnx, iys:iye), ty(nnx, iys:iye, izs:ize)
    REAL :: flux_u(nnx, iys:iye), flux_l(nnx, iys:iye)
    REAL :: taut3_u(nnx, iys:iye, nscl), taut3_l(nnx, iys:iye, nscl)
    REAL :: Sc, tscal, kconst

    ! SET SIGN FOR OCEAN SIMULATIONS THAT USE MONOTONE
    sgn = -1.0
    upwn = 2.0
    IF (iupwnd /= 1) upwn = 1.0

    ! OUTER LOOP OVER Z
    DO iz = izs, ize
        izm2 = iz - 2
        izm1 = iz - 1
        izp1 = iz + 1
        izp2 = iz + 2
        weit = dzw(iz) / (dzw(iz) + dzw(izp1))
        weit1 = 1.0 - weit
        weit3 = dzw(izm1) / (dzw(iz) + dzw(izm1))
        weit4 = 1.0 - weit3
        dzw2_i = 1.0 / (dzw(iz) + dzw(izp1))
        dzw3_i = 2.0 * dzw2_i
        DO iy = iys, iye
            DO ix = 1, nnx
                tx(ix, iy) = t(ix, iy, iscl, iz)
            END DO
        END DO

        CALL xderivp(tx(1, iys), trigx(1, 1), xk(1), nnx, iys, iye)

        ! COMPUTE TAU_T3 AT IZ-1
        IF (iz /= 1 .OR. ibcl /= 0) THEN
            DO iy = iys, iye
                DO ix = 1, nnx
                    taut3_l(ix, iy, iscl) = -vis_sv(ix, iy, izm1) * (t(ix, iy, iscl, iz) - &
                                                                     t(ix, iy, iscl, izm1)) * dzu_i(iz)
                END DO
            END DO
        ELSE
            DO iy = iys, iye
                DO ix = 1, nnx
                    taut3_l(ix, iy, iscl) = taut3m(ix, iy, iscl)
                END DO
            END DO
        END IF

        ! SGS TAU_T1,_T3 AND RESOLVED U*THETA SCALAR FLUXES
        ! SKEW SYMMETRIC ADVECTIVE TERM 0.5(UDT/DX + D/DX(UT))
        DO iy = iys, iye
            DO ix = 1, nnx
                taut3_u(ix, iy, iscl) = -vis_sv(ix, iy, iz) * (t(ix, iy, iscl, izp1) - &
                                                               t(ix, iy, iscl, iz)) * dzu_i(izp1)
                fnt1(ix, iy, iz) = -0.5 * ((vis_s(ix, iy, iz) + vis_s(ix, iy, izm1)) * &
                                           tx(ix, iy) - upwn * t(ix, iy, iscl, iz) * (u(ix, iy, iz) + stokes(iz) * &
                                                                                      dir_x))
            END DO
        END DO

        CALL xderivp(fnt1(1, iys, iz), trigx(1, 1), xk(1), nnx, iys, iye)

        DO iy = iys, iye
            DO ix = 1, nnx
                r4(ix, iy, iscl, iz) = -fnt1(ix, iy, iz) - (taut3_u(ix, iy, iscl) - &
                                                            taut3_l(ix, iy, iscl)) * dzw_i(iz)
            END DO
        END DO

        IF (iupwnd /= 1) THEN

            ! SKEW SYMMETRIC ADVECTIVE FORM FOR VERTICAL FLUX = 0.5(WDT/DZ + D/DZ(WT))
            DO iy = iys, iye
                DO ix = 1, nnx
                    theta_u = weit1 * t(ix, iy, iscl, iz) + weit * t(ix, iy, iscl, izp1)
                    theta_l = weit3 * t(ix, iy, iscl, iz) + weit4 * t(ix, iy, iscl, izm1)
                    r4(ix, iy, iscl, iz) = r4(ix, iy, iscl, iz) - 0.5 * (u(ix, iy, iz) + &
                                                                         stokes(iz) * dir_x) * tx(ix, iy) - 0.5 * (w(ix, iy, iz) * theta_u - &
                                                                                                                   w(ix, iy, izm1) * theta_l) * dzw_i(iz)
                    r4(ix, iy, iscl, iz) = r4(ix, iy, iscl, iz) - 0.25 * (w(ix, iy, iz) * &
                                                                          (t(ix, iy, iscl, izp1) - t(ix, iy, iscl, iz)) * dzu_i(izp1) + &
                                                                          w(ix, iy, izm1) * (t(ix, iy, iscl, iz) - t(ix, iy, iscl, izm1)) * &
                                                                          dzu_i(iz))
                END DO
            END DO
        ELSE
            ! Z-DIRECTION SPECIAL
            IF (iz == 1) THEN
                DO iy = iys, iye
                    DO ix = 1, nnx

                        ! AIR-SEA FLUX BC
                        IF (iscl == 2) THEN
                            Sc = 0.0d0
                            kconst = 0.0d0
                            kbub = 0.0d0
                            tscal = 0.0d0
                            tscal = t(ix, iy, 1, iz) - 273.15d0

                            ! CALCULATE SCHMIDT NUMBER (WANNINKOF, 1992)
                            Sc = 2073.1d0 - 125.62d0 * tscal + 3.6276d0 * tscal * tscal - &
                                 0.043219d0 * tscal * tscal * tscal

                            ! CALCULATE PISTON VELOCITY (WANNINKOF, 1992)
                            kconst = (2.77778d-6) * 0.31d0 * ws10 * ws10 * SQRT(660.0d0 / Sc)

                            ! CALCULATE BUBBLE PARAMETERIZATION
                            fug = 400.0 * exp(0.0413 * (tscal - 20.0)) * 1.0e-6   ! (WANNINKHOF, 2022)
                            khen = t(ix, iy, iscl, iz) * (10**(-6)) / fug         ! (EMERSON & HAMME, 2022)
                            bet_ost = khen * Rgas * t(ix, iy, 1, iz) * (1 / 0.101325)! (EMERSON & HAMME, 2022)
                            u_tau = SQRT(1.0 * (8.5e-4) * ws10 * ws10 / 1000.0)
                            wh = 0.0246 * ws10 * ws10; 
                            wa = 4.02e-7 * (u_tau * wh / (1.46e-5))**0.96
                            kbub = (2450.0 * wa) / (bet_ost * (1.0 + (14.0 * bet_ost * Sc**(-0.5))** &
                                                               (-1.0 / 1.2))**1.2)                       ! (WOOLF, 1997)

                            ! CALCULATE SURFACE FLUX RATE
                            IF (co2_asflux == 1) THEN
                                flux_l(ix, iy) = (kconst) * ((c1 * 1.10) - t(ix, iy, iscl, iz))
                            ELSEIF (co2_asflux == 2) THEN
                                flux_l(ix, iy) = (kconst + kbub) * ((c1 * 1.10) - t(ix, iy, iscl, iz))
                            ELSE
                                flux_l(ix, iy) = 0
                            END IF

                        ELSE
                            flux_l(ix, iy) = sgn * 0.5 * w(ix, iy, izm1) * (t(ix, iy, iscl, izm1) + &
                                                                            t(ix, iy, iscl, iz))
                        END IF

                        flux_u(ix, iy) = AMAX1(sgn * w(ix, iy, iz), 0.) * (t(ix, iy, iscl, iz) + &
                                                                           rlim(t(ix, iy, iscl, izp1), t(ix, iy, iscl, iz), &
                                                                                t(ix, iy, iscl, izm1))) + AMIN1(sgn * w(ix, iy, iz), 0.) * &
                                         (t(ix, iy, iscl, izp1) + rlim(t(ix, iy, iscl, iz), &
                                                                       t(ix, iy, iscl, izp1), t(ix, iy, iscl, izp2)))
                    END DO
                END DO
            ELSE IF (iz == nnz) THEN
                DO iy = iys, iye
                    DO ix = 1, nnx
                        flux_u(ix, iy) = sgn * 0.5 * w(ix, iy, iz) * (t(ix, iy, iscl, izp1) + &
                                                                      t(ix, iy, iscl, iz))
                        flux_l(ix, iy) = AMAX1(sgn * w(ix, iy, izm1), 0.) * &
                                         (t(ix, iy, iscl, izm1) + rlim(t(ix, iy, iscl, iz), &
                                                                       t(ix, iy, iscl, izm1), t(ix, iy, iscl, izm2))) + AMIN1(sgn * &
                                                                                                                              w(ix, iy, izm1), 0.) * (t(ix, iy, iscl, iz) + rlim( &
                                                                                                                                                      t(ix, iy, iscl, izm1), t(ix, iy, iscl, iz), t(ix, iy, iscl, izp1)))
                    END DO
                END DO
            ELSE
                DO iy = iys, iye
                    DO ix = 1, nnx
                        flux_u(ix, iy) = AMAX1(sgn * w(ix, iy, iz), 0.) * (t(ix, iy, iscl, iz) + &
                                                                           rlim(t(ix, iy, iscl, izp1), t(ix, iy, iscl, iz), &
                                                                                t(ix, iy, iscl, izm1))) + AMIN1(sgn * w(ix, iy, iz), 0.) * &
                                         (t(ix, iy, iscl, izp1) + rlim(t(ix, iy, iscl, iz), &
                                                                       t(ix, iy, iscl, izp1), t(ix, iy, iscl, izp2)))
                        flux_l(ix, iy) = AMAX1(sgn * w(ix, iy, izm1), 0.) * &
                                         (t(ix, iy, iscl, izm1) + rlim(t(ix, iy, iscl, iz), &
                                                                       t(ix, iy, iscl, izm1), t(ix, iy, iscl, izm2))) + AMIN1(sgn * &
                                                                                                                              w(ix, iy, izm1), 0.) * (t(ix, iy, iscl, iz) + rlim( &
                                                                                                                                                      t(ix, iy, iscl, izm1), t(ix, iy, iscl, iz), t(ix, iy, iscl, izp1)))
                    END DO
                END DO
            END IF

            ! SUM VERTICAL MONOTONE FLUX
            DO iy = iys, iye
                DO ix = 1, nnx
                    r4(ix, iy, iscl, iz) = r4(ix, iy, iscl, iz) - sgn * (flux_u(ix, iy) - &
                                                                         flux_l(ix, iy)) * dzw_i(iz)
                END DO
            END DO
        END IF

        ! SAVE SGS FLUXES FOR PRINTOUT, GATHER SUMS ON EXIT
        IF (istep == 1) THEN
            utsb(iz, iscl) = 0.0
            wtsb(iz, iscl) = 0.0
            DO iy = iys, iye
                DO ix = 1, nnx
                    wtsb(iz, iscl) = wtsb(iz, iscl) + taut3_u(ix, iy, iscl)
                    utsb(iz, iscl) = utsb(iz, iscl) - 0.5 * (vis_s(ix, iy, iz) + &
                                                             vis_s(ix, iy, izm1)) * tx(ix, iy)
                END DO
            END DO
            utsb(iz, iscl) = utsb(iz, iscl) * fnxy
            wtsb(iz, iscl) = wtsb(iz, iscl) * fnxy
        END IF
    END DO

    ! OUTER LOOP OVER Z FOR Y-DEPENDENCE
    DO iz = izs, ize
        DO iy = iys, iye
            DO ix = 1, nnx
                ty(ix, iy, iz) = t(ix, iy, iscl, iz)
            END DO
        END DO
    END DO

    ! Y DERIVATIVE OF T FOR [IZS:IZE]
    CALL yd_mpi(ty(1, iys, izs), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs, ize, myid, ncpu_s, numprocs)

    ! ADD SKEW SYMMETRIC ADVECTIVE FLUX AND SGS FLUX TO Y-DERIV COMPUTATION.
    ! CHECK FOR MONOTONE
    DO iz = izs, ize
        izm1 = iz - 1
        DO iy = iys, iye
            DO ix = 1, nnx
                fnt1(ix, iy, iz) = -0.5 * ((vis_s(ix, iy, iz) + vis_s(ix, iy, izm1)) * &
                                           ty(ix, iy, iz) - upwn * t(ix, iy, iscl, iz) * (v(ix, iy, iz) + &
                                                                                          stokes(iz) * dir_y))
            END DO
        END DO

        IF (iupwnd /= 1) THEN
            DO iy = iys, iye
                DO ix = 1, nnx
                    r4(ix, iy, iscl, iz) = r4(ix, iy, iscl, iz) - 0.5 * (v(ix, iy, iz) + &
                                                                         stokes(iz) * dir_y) * ty(ix, iy, iz)
                END DO
            END DO
        END IF
    END DO

    ! Y- DERIVATIVES OF SCALAR FLUXES FOR [IZS:IZE]
    CALL yd_mpi(fnt1(1, iys, izs), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs, ize, myid, ncpu_s, numprocs)

    DO iz = izs, ize
        DO iy = iys, iye
            DO ix = 1, nnx
                r4(ix, iy, iscl, iz) = r4(ix, iy, iscl, iz) - fnt1(ix, iy, iz)
            END DO
        END DO
    END DO

    ! SAVE SGS FLUXES FOR PRINTOUT
    IF (istep == 1) THEN
        DO iz = izs, ize
            vtsb(iz, iscl) = 0.0
            DO iy = iys, iye
                DO ix = 1, nnx
                    vtsb(iz, iscl) = vtsb(iz, iscl) - 0.5 * (vis_s(ix, iy, iz) + &
                                                             vis_s(ix, iy, izm1)) * ty(ix, iy, iz)
                END DO
            END DO
            vtsb(iz, iscl) = vtsb(iz, iscl) * fnxy
        END DO
    END IF

    RETURN
END SUBROUTINE
