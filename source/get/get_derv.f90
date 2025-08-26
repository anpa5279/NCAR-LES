SUBROUTINE get_derv
! GET UX,UY,VX,VY AT ALL Z FOR THIS NODE USING PARALLEL FFT
! CAN BE IMPROVED USING EXCHANGE TO SEND DERIVATIVES

    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats

    iz_ss = izs - 1
    iz_ee = ize + 1
    IF (iss == 0) THEN
        iz_ss = izs
    END IF

    IF (ise == numprocs - 1) THEN
        iz_ee = ize
    END IF

    ! MAKE SURE <W> = 0
    DO iz = izs - 1, ize + 1
        w_sum = 0.0
        DO iy = iys, iye
            DO ix = 1, nnx
                w_sum = w_sum + w(ix, iy, iz)
            END DO
        END DO
        w_sum = w_sum * fnxy

        CALL mpi_sum_xy(w_sum, myid, iss, ise, 1)

        DO iy = iys, iye
            DO ix = 1, nnx
                w(ix, iy, iz) = w(ix, iy, iz) - w_sum
            END DO
        END DO
    END DO

    DO iz = izs - 1, ize + 1
        DO iy = iys, iye
            DO ix = 1, nnx
                ux(ix, iy, iz) = u(ix, iy, iz)
                vx(ix, iy, iz) = v(ix, iy, iz)
                wx(ix, iy, iz) = w(ix, iy, iz)
                uy(ix, iy, iz) = u(ix, iy, iz)
                vy(ix, iy, iz) = v(ix, iy, iz)
                wy(ix, iy, iz) = w(ix, iy, iz)
            END DO
        END DO

        CALL xderivp(ux(1, iys, iz), trigx(1, 1), xk(1), nnx, iys, iye)
        CALL xderivp(vx(1, iys, iz), trigx(1, 1), xk(1), nnx, iys, iye)
        CALL xderivp(wx(1, iys, iz), trigx(1, 1), xk(1), nnx, iys, iye)
    END DO

    ! GET Y DERIVATIVES FOR (U,V,W)
    CALL yd_mpi(uy(1, iys, izs - 1), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs - 1, ize + 1, myid, ncpu_s, numprocs)
    CALL yd_mpi(vy(1, iys, izs - 1), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs - 1, ize + 1, myid, ncpu_s, numprocs)
    CALL yd_mpi(wy(1, iys, izs - 1), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs - 1, ize + 1, myid, ncpu_s, numprocs)

    RETURN
END SUBROUTINE
