SUBROUTINE comp2
! ADD P GRADIENTS TO RHS
! USE ALREADY DEFINED P AT IZE+1 TO GET W (SEE SR. PRESSURE)

    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats

    REAL :: fnt1(nnx, iys:iye, izs:ize), fnt2(nnx, iys:iye)
    REAL :: r3_sum(1:nnz)

    INCLUDE 'mpif.h'

    INTEGER :: istatus(mpi_status_size)

    DO iz = 1, nnz
        r3_sum(iz) = 0.0
    END DO

    ! DP/DY AT ALL Z
    DO iz = izs, ize
        DO iy = iys, iye
            DO ix = 1, nnx
                fnt1(ix, iy, iz) = p(ix, iy, iz)
            END DO
        END DO
    END DO

    CALL yd_mpi(fnt1(1, iys, izs), trigx(1, 2), yk(1), nnx, nny, ixs, ixe, ix_s, ix_e, &
                iys, iye, iy_s, iy_e, izs, ize, myid, ncpu_s, numprocs)

    DO iz = izs, ize
        izm1 = iz - 1
        izp1 = iz + 1
        DO iy = iys, iye
            DO ix = 1, nnx
                fnt2(ix, iy) = p(ix, iy, iz)
            END DO
        END DO

        CALL xderivp(fnt2(1, iys), trigx(1, 1), xk(1), nnx, iys, iye)

        DO iy = iys, iye
            DO ix = 1, nnx
                r1(ix, iy, iz) = r1(ix, iy, iz) - fnt2(ix, iy)
                r2(ix, iy, iz) = r2(ix, iy, iz) - fnt1(ix, iy, iz)
            END DO
        END DO

        IF (iz /= nnz) THEN
            DO iy = iys, iye
                DO ix = 1, nnx
                    r3(ix, iy, iz) = r3(ix, iy, iz) - (p(ix, iy, izp1) - p(ix, iy, iz)) * dzu_i(izp1)
                    r3_sum(iz) = r3_sum(iz) + r3(ix, iy, iz)
                END DO
            END DO
            r3_sum(iz) = r3_sum(iz) * fnxy
        END IF

        ! TIME STEPPING WITH 3RD ORDER RK METHOD FIRST W VARIABLES
        IF (iz /= nnz) THEN
            DO iy = iys, iye
                DO ix = 1, nnx
                    e(ix, iy, iz) = e(ix, iy, iz) + dtgama * r5(ix, iy, iz)
                END DO
            END DO
        ELSE

            ! UPDATE WOUT AND EOUT BY SETTING = TO BC VALUES
            DO iy = iys, iye
                DO ix = 1, nnx
                    w(ix, iy, iz) = wbc(ix, iy, 1)
                    e(ix, iy, iz) = ebc(ix, iy, 1)
                    r3(ix, iy, iz) = 0.0
                    r5(ix, iy, iz) = 0.0
                END DO
            END DO
        END IF

        ! NOW ALL U-VARIABLES
        DO iy = iys, iye
            DO ix = 1, nnx
                u(ix, iy, iz) = u(ix, iy, iz) + dtgama * r1(ix, iy, iz)
                v(ix, iy, iz) = v(ix, iy, iz) + dtgama * r2(ix, iy, iz)
            END DO
        END DO

        DO iscl = 1, nscl
            DO iy = iys, iye
                DO ix = 1, nnx
                    t(ix, iy, iscl, iz) = t(ix, iy, iscl, iz) + dtgama * r4(ix, iy, iscl, iz)
                END DO
            END DO
        END DO
    END DO

    ! GATHER PARTIAL SUMS FOR W COMPUTATION
    CALL mpi_sum_z(r3_sum, i_root, myid, nnz, 1)

    DO iz = izs, MIN(ize, nnz - 1)
        DO iy = iys, iye
            DO ix = 1, nnx
                r3(ix, iy, iz) = r3(ix, iy, iz) - r3_sum(iz)
                w(ix, iy, iz) = w(ix, iy, iz) + dtgama * r3(ix, iy, iz)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
