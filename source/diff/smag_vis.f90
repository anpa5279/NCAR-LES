SUBROUTINE smag_vis(istep)
    !! Classic Smagorinsky LES model.

    USE pars
    USE fields
    USE con_data
    USE con_stats

    IMPLICIT NONE

    REAL :: wz, uz, vz, s11, s22, s33, s12, s13, s23, sij2, delta
    REAL :: d_im, d_ip, d_i, wt_l, wt_h
    INTEGER :: istep, iz, i, j, izm1, izp1

    ! Step 0: Change the Lower BC for u and v
    ! -------
    ! We want zero vertical SGS flux across the boundaries, which requires that
    ! uz = vz = 0.0 at the z-face pt, which affects the gradient computation at
    ! the first interior z-center pt for u and v. Therefore override the BCs set
    ! by Lower.f90 for right now.
    ! I believe Lower.f90 sets u = v = 0 for an atmospheric BL still,
    ! instead of uz = vz = 0 for an oceanic BL, which no one caught because it
    ! didn't matter when using a staggered vertical grid and one-equation SGS
    ! model. I could be wrong though, and maybe that BC is correct, and I need
    ! to rethink this...
    IF (izs == 1) THEN
        u(:, :, 0) = u(:, :, 1)
        v(:, :, 0) = v(:, :, 1)
    END IF

    ! Step 1: Compute SGS viscosity at z-center interior points
    ! -------
    ! NOTE: z-center points are always halfway between z-face points, therefore, to
    ! interpolate from face to center, it's a simple (equal-weighted) average.
    DO iz = izs, ize
        izm1 = iz - 1
        izp1 = iz + 1
        delta = dsl_z(iz)

        ! 2nd-order finite difference coefficients for z-center gradient from
        ! z-center data on a non-uniform grid
        d_im = 1.0 / (dzw(iz) + dzw(izm1)) ! = 1 / 2dz on uniform grid
        d_ip = 1.0 / (dzw(izp1) + dzw(iz)) ! = 1 / 2dz on uniform grid
        d_i = d_im - d_ip                  ! = 0 on uniform grid

        DO j = iys, iye
            DO i = 1, nnx
                ! Calculate vertical gradients at z-center points
                wz = (w(i, j, iz) - w(i, j, izm1)) * dzw_i(iz)
                uz = d_ip * u(1, j, izp1) + d_i * u(1, j, iz) - d_im * u(1, j, izm1)
                vz = d_ip * v(1, j, izp1) + d_i * v(1, j, iz) - d_im * v(1, j, izm1)

                ! Strain rate tensor terms squared, sij = [1/2*(u_<i,j> + u_<j,i>)]
                s11 = ux(i, j, iz)**2
                s22 = vy(i, j, iz)**2
                s33 = wz**2
                s12 = 0.25 * (uy(i, j, iz) + vx(i, j, iz))**2
                s13 = 0.25 * (uz + 0.5 * (wx(i, j, iz) + wx(i, j, izm1)))**2
                s23 = 0.25 * (vz + 0.5 * (wy(i, j, iz) + wy(i, j, izm1)))**2
                sij2 = s11 + s22 + s33 + 2.0 * (s12 + s13 + s23)

                ! using vis_sv to store z-center viscosity values
                vis_sv(i, j, iz) = (csmag * delta)**2 * SQRT(2.0 * sij2)
                IF (isnan(vis_sv(i, j, iz))) THEN
                    IF (l_root) WRITE (6, *) "wt_h=", wt_h, "wt_l=", wt_l
                    IF (l_root) WRITE (6, *) "dzw(iz)=", dzw(iz), "dzw(izp1)=", dzw(izp1), "dzw_i(iz) = ", dzw_i(iz), "dzu_i(iz)=", dzu_i(iz)
                    IF (l_root) WRITE (6, *) "dzw=", dzw
                    STOP
                END IF
            END DO
        END DO
    END DO

    ! Step 2: Interpolate to z-face points, and add appropriate BCs
    ! -------
    call smag_exchange ! fills halo cells just for vis_sv

    ! Step 3: Interpolate to z-face points
    ! -------
    DO iz = izs, ize
        izp1 = iz + 1
        wt_h = dzw(iz) / (dzw(iz) + dzw(izp1))
        wt_l = 1.0 - wt_h

        DO j = iys, iye
            DO i = 1, nnx
                vis_m(i, j, iz) = wt_l * vis_sv(i, j, iz) + wt_h * vis_sv(i, j, izp1)
            END DO
        END DO
    END DO

    ! Step 4: Set the Boundary Conditions
    ! -------
    ! Lower BC (air-sea interface):
    ! No SGS diffusion of quantities across interface
    ! -- RHS_X override tau13, tau23, and taut3 with only the flux BC
    !    and ignore SGS stresses, which covers this requirement.
    !
    ! Horizontal SGS diffusion at first internal point should depend only on the
    ! viscosity at that internal point.
    ! -- This requires that when RHS_X averages the face-located vis_x, the
    !    result is exactly what was computed in Step 1 for the center-located
    !    vis_x
    !    vis_c(1) = 0.5 * (vis_f(0) + vis_f(1))
    !    -> vis_f(0) = 2.0 * vis_c(1) - vis_f(1)
    ! -- This is different than TKE-based model, which sets vis_f(0) = vis_f(1),
    !    resulting in vis_c(1) = vis_f(1), which isn't perfectly accurate.
    IF (izs == 1) THEN
        vis_m(:, :, izs-1) = 2.0 * vis_sv(:, :, izs) - vis_m(:, :, izs)
    END IF

    ! The Upper BC isn't exactly a reflection/slip wall BC.
    ! Temperature (and scalars) can have a non-zero gradient. The TKE-based
    ! SGS model relies on e(nnz) = 0.0 to zero out all SGS diffusion
    ! at the boundary, which distorts horizontal SGS diffusion at the last
    ! interior point for u, v, and t. Here we have the option to better match
    ! velocities and scalars separately than in the TKE-based model.
    !
    ! No vertical SGS diffusion
    ! -- tau13 = tau23 = 0 at face-pt nnz because all velocity gradients are
    !    already 0. But, dT/dz is not 0, so we need to make vis_sv = 0.
    !
    ! Properly-matched horizontal SGS diffusion
    ! -- same procedure as above for vis_m and vis_s
    IF (ize == nnz) THEN
        vis_m(:, :, ize) = 2.0 * vis_sv(:, :, ize) - vis_m(:, :, ize-1)
        vis_s(:, :, ize) = 3.0 * vis_m(:, :, ize)
        vis_sv(:, :, ize) = 0.0
    ELSE
        vis_s(:, :, ize) = 3.0 * vis_m(:, :, ize)
        vis_sv(:, :, ize) = vis_s(:, :, ize)
    END IF

    ! Step 5: Update the rest of vis_s and vis_sv
    ! -------
    vis_s(:, :, izs-1:ize-1) = 3.0 * vis_m(:, :, izs-1:ize-1)
    vis_sv(:, :, izs-1:ize-1) = vis_s(:, :, izs-1:ize-1)

    ! Step 6: Revert the change to the Lower BC for u and v
    ! -------
    ! I'm not willing to bet my life that there isn't something somewhere in the
    ! code that needs the u = v = 0 BC to be correct, which is why we're keeping
    ! Lower.f90 as-is and doing this thing.
    IF (izs == 1) THEN
        u(:, :, 0) = ubc(:, :, 2)
        v(:, :, 0) = vbc(:, :, 2)
    END IF

    RETURN
END SUBROUTINE
