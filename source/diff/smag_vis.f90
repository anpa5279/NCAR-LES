SUBROUTINE smag_vis(istep)
! Classic Smagorinsky LES model.
! WARNING: assumes ivis0 == 0 (set in pars.f90).

    USE pars
    USE fields
    USE con_data
    USE con_stats

    REAL :: d_grid, wz, uz, vz, wzp, sij2, wt_h, wt_l
    INTEGER :: iz, i, j, izm1, izp1

    DO iz = izs - 1, ize + 1
        izm1 = iz - 1
        izp1 = iz + 1
        wt_h = dzw(iz) / (dzw(iz) + dzw(izp1))
        wt_l = 1.0 - wt_h
        d_grid = dsl_z(iz) ! dsl_z gives correct SGS length scale given dealiased x and y derivatives

        DO j = iys, iye
            DO i = 1, nnx
                ! Calculate vertical gradients.
                wz = (w(i, j, iz) - w(i, j, izm1)) * dzw_i(iz)
                uz = (u(i, j, iz) - u(i, j, izm1)) * dzu_i(iz)
                vz = (v(i, j, iz) - v(i, j, izm1)) * dzu_i(iz)
                wzp = (w(i, j, izp1) - w(i, j, iz)) * dzw_i(izp1)

                ! Strain rate tensor terms squared. strain rate tens: sij= 1/2*(u_<i,j> + u_<j,i>)
                s11 = wt_l * ux(i, j, iz)**2 + wt_h * ux(i, j, izp1)**2
                s22 = wt_l * vy(i, j, iz)**2 + wt_h * vy(i, j, izp1)**2
                s33 = wt_h * wzp**2 + wt_l * wz**2
                s12 = wt_l * (0.5 * (uy(i, j, iz) + vx(i, j, iz)))**2 &
                      + wt_h * (0.5 * (uy(i, j, izp1) + vx(i, j, izp1)))**2

                ! Note from Colin: why s13 and s23 are not interpolated between iz and izp1
                !   seems to be related to differing cell center/face alignments in RHS_UVW.
                !   Either way, let's make sure we match TKE_VIS here (less the xy-means).
                s13 = (0.5 * (uz + wx(i, j, iz)))**2
                s23 = (0.5 * (vz + wy(i, j, iz)))**2

                sij2 = s11 + s22 + s33 + 2.0 * (s12 + s13 + s23)
                vis_m(i, j, iz) = (csmag * d_grid)**2 * SQRT(2.0 * sij2)
                vis_s(i, j, iz) = 3.0 * vis_m(i, j, iz) 
                vis_sv(i, j, iz) = vis_s(i, j, iz)
                IF (isnan(vis_m(i, j, iz))) THEN ! debugging, can delete later
                    IF (l_root) WRITE (6, *) "NaNs appeared in smag_vis at i=",i, " j=",j, " iz=",iz
                    IF (l_root) WRITE (6, *) "wt_h=", wt_h, "wt_l=", wt_l
                    IF (l_root) WRITE (6, *) "dzw(iz)=", dzw(iz), "dzw(izp1)=", dzw(izp1), "dzw_i(iz) = ", dzw_i(iz), "dzu_i(iz)=", dzu_i(iz)
                    IF (l_root) WRITE (6, *) "dzw=", dzw
                    CALL mpi_finalize(ierr)
                    STOP
                END IF
            END DO
        END DO
    END DO
    RETURN
END SUBROUTINE
