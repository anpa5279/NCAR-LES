SUBROUTINE get_max
! ROUTINE COMPUTES MAX VELOCITIES AS SWEEP THROUGH THE VELOCITY FIELD

    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats
    INCLUDE 'mpif.h'

    REAL :: u_send(6), u_recv(6)

    dx_i = 1.0 / dx
    dy_i = 1.0 / dy

    u_temp = 0.0
    v_temp = 0.0
    w_temp = 0.0
    e_temp = 0.0
    vis_temp = 0.0

    DO iz = izs, ize
        u_xy = 0.0
        v_xy = 0.0
        w_xy = 0.0
        e_xy = 0.0
        DO iy = iys, iye
            DO ix = 1, nnx
                u_xy = AMAX1(u_xy, ABS(u(ix, iy, iz) + stokes(iz) * dir_x))
                v_xy = AMAX1(v_xy, ABS(v(ix, iy, iz) + stokes(iz) * dir_y))
                w_xy = AMAX1(w_xy, ABS(w(ix, iy, iz)))
                e_xy = AMAX1(e_xy, e(ix, iy, iz))
            END DO
        END DO

        u_xy = u_xy * dx_i
        v_xy = v_xy * dy_i
        wsav = w_xy
        w_xy = w_xy / ABS(dzw(iz))
        vis_xy = 3.0 * ck * dsl_z(iz) * SQRT(e_xy) / AMIN1(dx, dy, dzw(iz))**2

        u_temp = AMAX1(u_xy, u_temp)
        v_temp = AMAX1(v_xy, v_temp)
        w_temp = AMAX1(w_xy, w_temp)
        e_temp = AMAX1(e_xy, e_temp)
        vis_temp = AMAX1(vis_xy, vis_temp)
    END DO

    u_send(1) = u_temp
    u_send(2) = v_temp
    u_send(3) = w_temp
    u_send(4) = wsav
    u_send(5) = e_temp
    u_send(6) = vis_temp

    CALL mpi_allreduce(u_send, u_recv, 6, mpi_real8, mpi_max, mpi_comm_world, ierror)

    umax = u_recv(1)
    vmax = u_recv(2)
    wmax = u_recv(3)
    wABS = u_recv(4)
    emax = u_recv(5)
    vismax = u_recv(6)

    RETURN
END SUBROUTINE
