SUBROUTINE smag_vis(istep)
! Classic Smagorinsky LES model.
! WARNING: assumes ivis0 == 0 (set in pars.f90).

    USE pars
    USE fields
    USE con_data
    USE con_stats

    IMPLICIT NONE

    REAL :: wz, uzp, vzp, wzp, s11, s22, s33, s12, s13, s23, sij2, weit, weit1, d_grid
    INTEGER :: istep, iz, i, j, izm1, izp1

    ! Index (i, j, izs-2) is out of bounds when computing wz for izs-1, and (i, j, ize+2) is out
    ! out bounds when computing uzp, vzp, wzp for ize+1. Therefore, we must supply vis_x values
    ! for izs-1 and ize+1 separately, regardless of whether izs=1 or ize=nnz (BCs). This is simply
    ! a halo-cell inadequacy in NCARLES for computing Smagorinsky viscosity, and it implies that
    ! the solution will change as the MPI decomposition is changed (i.e., solution becomes
    ! dependent on values of ncpu_s and numprocs), which is not good.
    DO iz = izs, ize
        izm1 = iz - 1
        izp1 = iz + 1
        weit = dzw(iz) / (dzw(iz) + dzw(izp1))
        weit1 = 1.0 - weit
        d_grid = dsl_z(iz) ! dsl_z gives correct SGS length scale given dealiased x and y derivatives

        DO j = iys, iye
            DO i = 1, nnx
                ! Calculate vertical gradients.
                wz = (w(i, j, iz) - w(i, j, izm1)) * dzw_i(iz)
                uzp = (u(i, j, izp1) - u(i, j, iz)) * dzu_i(izp1)
                vzp = (v(i, j, izp1) - v(i, j, iz)) * dzu_i(izp1)
                wzp = (w(i, j, izp1) - w(i, j, iz)) * dzw_i(izp1)

                ! Strain rate tensor terms squared. strain rate tens: sij= 1/2*(u_<i,j> + u_<j,i>)
                s11 = weit1 * ux(i, j, iz)**2 + weit * ux(i, j, izp1)**2
                s22 = weit1 * vy(i, j, iz)**2 + weit * vy(i, j, izp1)**2
                s33 = weit * wzp**2 + weit1 * wz**2
                s12 = weit1 * (0.5 * (uy(i, j, iz) + vx(i, j, iz)))**2 &
                      + weit * (0.5 * (uy(i, j, izp1) + vx(i, j, izp1)))**2

                s13 = (0.5 * (uzp + wx(i, j, iz)))**2
                s23 = (0.5 * (vzp + wy(i, j, iz)))**2

                sij2 = s11 + s22 + s33 + 2.0 * (s12 + s13 + s23)
                vis_m(i, j, iz) = (csmag * d_grid)**2 * SQRT(2.0 * sij2)
                vis_s(i, j, iz) = 3.0 * vis_m(i, j, iz) ! 1 + 2 * (alk/dslk) = 3 for non-stability-corrected length scales
                vis_sv(i, j, iz) = vis_s(i, j, iz)
            END DO
        END DO
    END DO

    ! Lower BC sets these stored values @iz=0:
    !   r3 = 0.0
    !   ux = uy = vx = vy = wx = wy = 0
    !   tau_13 = wind stress from SUFTO
    !   tau_23 = wind stress from SUFTO
    ! Therefore, tau_11 = tau_22 = tau_12 = 0 identically, tau_33 is ignored (r3 = 0), and tau_13
    ! and tau_23 do not require vis_m. However, the vis_x @iz=0 are added to vis_x @iz=1 to
    ! compute tau_11, tau_12, tau_22, and the tau_t's @iz=1. Therefore the best thing to do is set
    ! them equal to the vis_x @iz=1, as in the TKE_VIS routines.
    IF (izs == 1) THEN
        DO j = iys, iye
            DO i = 1, nnx
                vis_m(i, j, izs - 1) = vis_m(i, j, izs)
                vis_s(i, j, izs - 1) = vis_s(i, j, izs)
                vis_sv(i, j, izs - 1) = vis_sv(i, j, izs)
            END DO
        END DO
    END IF

    ! Upper BC sets these stored values @iz=nnz+1:
    !   u = u(nnz) => implies `uzp` @(nnz:nnz+1) = 0
    !   v = v(nnz) => implies `vzp` @(nnz:nnz+1) = 0
    !   ux = uy = vx = vy = wx = wy = 0
    ! and @iz=(nnz:nnz+1) sets these values:
    !   w = 0.0 => implies `wzp` @(nnz) = 0, `wz` @(nnz+1) = 0
    !   r3 = 0.0.
    ! Therefore, all sij = 0.0 @nnz+1 identically.
    ! This means the tau's and vis_x's should all be 0.0 @ nnz+1. In fact, the
    ! vis_s/sv(ize+1) index is not even used in RHS_SCL, and thefore isn't needed.
    ! However, vis_m(ize+1) is only used in RHS_UVW
    ! where, if ize=nnz, the result is overwritten by r3 = 0.0. Therefore, we can just let whatever
    ! value got computed for ize+1 in the main DO loop above be left alone, as it does not
    ! computationally matter, just like in the TKE_VIS routines.
    IF (ize == nnz) THEN
        DO j = iys, iye
            DO i = 1, nnx
                vis_m(i, j, ize + 1) = 0.0
                vis_s(i, j, ize + 1) = 0.0
                vis_sv(i, j, ize + 1) = 0.0
            END DO
        END DO
    END IF

    ! SIDE NOTE: The unexpected Dirichlet BCs for u and v set in lower.f90 are irrelevant to
    ! rhs_uvw, as the gradients ux, uy, vx, vy, as well as the stresses tau_13 and tau_23, are all
    ! set separately. Therefore, only uz and vz are set by the u and v BC, and because of the
    ! staggered vertical grid, both of these gradients only appear advectively, multiplied by the
    ! w = 0.0 condition. My guess is this setup is a leftover from the old atmospheric physics.
    ! My expectation was that the "rigid lid" surface condition would be a free-slip wall, and
    ! therefore u and v would have gradient conditions (equaling wind-stress for u, zero for v).
    ! I don't know what ugal or ugcont are, but they are set to 0.0 in the code.

    RETURN
END SUBROUTINE
