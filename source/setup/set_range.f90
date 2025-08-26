SUBROUTINE set_range

    USE pars

    WRITE (nprt, 7002) nnx, nny, nnz

    ii = -1
    DO nn = 0, ncpu_z - 1
        CALL range(1, nnx + 2, ncpu_z, nn, lx_s, lx_e)
        CALL range(1, nnx, ncpu_z, nn, nx_s, nx_e)
        CALL range(1, nny, ncpu_z, nn, ly_s, ly_e)
        CALL range(1, nnz, ncpu_z, nn, mz_s, mz_e)

        DO mm = 0, ncpu_s - 1
            CALL range(1, nny, ncpu_s, mm, ny_s, ny_e)
            CALL range(1, nnx, ncpu_s, mm, nxy_s, nxy_e)
            CALL range(1, ncx, ncpu_s, mm, l2x_s, l2x_e)

            ii = ii + 1
            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1) * 2 + 1
            jx_e(ii) = l2x_e * 2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e

            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e

            iz_s(ii) = mz_s
            iz_e(ii) = mz_e

            is_s(ii) = (ii / ncpu_s) * ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
        END DO
    END DO

    iys = iy_s(myid)
    iye = iy_e(myid)
    jys = jy_s(myid)
    jye = jy_e(myid)
    ixs = ix_s(myid)
    ixe = ix_e(myid)
    jxs = jx_s(myid)
    jxe = jx_e(myid)
    kxs = kx_s(myid)
    kxe = kx_e(myid)
    mxs = mx_s(myid)
    mxe = mx_e(myid)
    izs = iz_s(myid)
    ize = iz_e(myid)

    ! GET STARTING AND ENDING PROCESSOR ID ON EACH VERTICAL SLAB
    iss = is_s(myid)
    ise = is_e(myid)

    ! DEBUG RANGES
    IF (l_debug) THEN
        WRITE (nprt, 1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn), jx_e(nn), &
                                  kx_s(nn), kx_e(nn), nn=0, numprocs - 1)
        WRITE (nprt, 1213) myid, (nn, iy_s(nn), iy_e(nn), jy_s(nn), jy_e(nn), &
                                  iz_s(nn), iz_e(nn), is_s(nn), is_e(nn), nn=0, numprocs - 1)
    END IF

    RETURN

! FORMAT
7002 FORMAT(' 7002 gridd nnx = ', i4, ' nny = ', i4, ' nnz = ', i4)
1200 FORMAT(' myid =  ', i4, /, ' nn', 5x, ' ixs ', 5x, ' ixe ', 5x, ' jxs ', 5x, &
           ' jxe ', 5x, ' kxs ', 5x, ' kxe', /, (7i6))
1213 FORMAT(' myid = ', i3, /, ' nn ', 3x, ' iys ', 5x, ' iye ', 5x, ' jys ', 5x, &
           ' jye ', 5x, ' izs ', 5x, ' ize', 5x, ' iss ', 5x, ' ise ', /, (9i6))

END SUBROUTINE
