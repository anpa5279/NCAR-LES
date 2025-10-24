SUBROUTINE smag_vis(istep)
! Classic Smagorinsky LES model.
! WARNING: assumes ivis0 == 0 (set in pars.f90).

    USE pars
    USE fields
    USE con_data
    USE con_stats

    REAL :: d_grid, wz, uz, vz, wzp, sij2
    INTEGER :: iz, i, j, izm1, izp1

    DO iz = izs - 1, ize + 1
        izm1 = iz - 1
        izp1 = iz + 1
        weit = dzw(iz) / (dzw(iz) + dzw(izp1))
        weit1 = 1.0 - weit
        d_grid = dsl_z(iz) ! dsl_z gives correct SGS length scale given dealiased x and y derivatives

        DO j = iys, iye
            DO i = 1, nnx
                ! Calculate vertical gradients.
                wz = (w(i, j, iz) - w(i, j, izm1)) * dzw_i(iz)
                uz = (u(i, j, iz) - u(i, j, izm1)) * dzu_i(iz)
                vz = (v(i, j, iz) - v(i, j, izm1)) * dzu_i(iz)
                wzp = (w(i, j, izp1) - w(i, j, iz)) * dzw_i(izp1)

                ! Strain rate tensor terms squared. strain rate tens: sij= 1/2*(u_<i,j> + u_<j,i>)
                s11 = weit1 * ux(i, j, iz)**2 + weit * ux(i, j, izp1)**2
                s22 = weit1 * vy(i, j, iz)**2 + weit * vy(i, j, izp1)**2
                s33 = weit * wzp**2 + weit1 * wz**2
                s12 = weit1 * (0.5 * (uy(i, j, iz) + vx(i, j, iz)))**2 &
                      + weit * (0.5 * (uy(i, j, izp1) + vx(i, j, izp1)))**2

                ! Note from Colin: why s13 and s23 are not interpolated between iz and izp1
                !   seems to be related to differing cell center/face alignments in RHS_UVW.
                !   Either way, let's make sure we match TKE_VIS here (less the xy-means).
                s13 = (0.5 * (uz + wx(i, j, iz)))**2
                s23 = (0.5 * (vz + wy(i, j, iz)))**2

                sij2 = s11 + s22 + s33 + 2.0 * (s12 + s13 + s23)
                vis_m(i, j, iz) = (csmag * d_grid)**2 * SQRT(2.0 * sij2)
                vis_s(i, j, iz) = 3.0 * vis_m(i, j, iz) 
                vis_sv(i, j, iz) = vis_s(i, j, iz)
                IF (isnan(vis_m(i, j, iz))) THEN
                    IF (l_root) WRITE (6, *) "NaNs appeared in smag_vis at i=",i, " j=",j, " iz=",iz
                    IF (l_root) WRITE (6, *) "weit=", weit, "weit1=", weit1
                    IF (l_root) WRITE (6, *) "dzw_i(iz) = ", dzw_i(iz), "dzu_i(iz)=", dzu_i(iz)
                    IF (l_root) WRITE (6, *) "uz=", uz, "vz=", vz, "wz=", wz, "wzp=", wzp, "ux(i, j, iz)=", ux(i, j, iz), "ux(i, j, izp1)=", ux(i, j, izp1), "vy(i, j, iz)=", vy(i, j, iz), "vy(i, j, izp1)=", vy(i, j, izp1)
                    IF (l_root) WRITE (6, *) "s11=",s11, "s22=", s22, "s33=", s33, "s23=", s23, "s13=", s13
                    IF (l_root) WRITE (6, *) "u(i, j, iz) = ", u(i, j, iz), "u(i, j, izm1)=", u(i, j, izm1), "v(i, j, iz)=", v(i, j, iz), "v(i, j, izm1)=", v(i, j, izm1), "w(i, j, iz)=", w(i, j, iz), "w(i, j, izm1)=", w(i, j, izm1), "w(i, j, izp1)=", w(i, j, izp1)
                    CALL mpi_finalize(ierr)
                    STOP
                END IF
            END DO
        END DO
    END DO
    RETURN
END SUBROUTINE
