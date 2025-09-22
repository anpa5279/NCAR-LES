SUBROUTINE vgrid(z_lo, z_ml, z_hi, nnz, z, l_root, ldebug)

    INTEGER :: nnz, iter
    REAL :: z(0:nnz + 1), z_lo, z_ml, z_hi, z_ratio, nnz_i, stretch, new_stretch, test
    LOGICAL :: l_root, l_debug

    real, parameter :: tol = 1.0e-10

    z_ratio = z_hi / z_lo
    nnz_i = 1.0 / FLOAT(nnz)
    stretch = 1.1   ! initial stretch factor to try
    iter = 0

    DO
        new_stretch = (z_ratio * (stretch - 1.0) + 1.0)**nnz_i
        test = ABS(1.0 - new_stretch / stretch)

        ! while loop exit condition
        if (test <= tol) exit

        ! while loop error condition
        iter = iter + 1
        IF (iter > 500) THEN
            IF (l_root) WRITE (6, 9000) stretch, new_stretch, iter
            STOP
        END IF
    END DO

    IF (l_root) WRITE (6, 9100) stretch, z_hi, z_lo, iter

    z(0) = 0.0
    z(1) = z_lo

    ! this will be NaN if stretch == 1.0
    DO iz = 2, nnz
        z(iz) = z_lo * (stretch**iz - 1.0) / (stretch - 1.0)
    END DO

    z(nnz) = z_hi
    z(nnz + 1) = 2.0 * z_hi - z(nnz - 1) ! z_hi + dz

    IF (l_root) WRITE (6, 5300) nnz

    RETURN

! FORMAT
9000 FORMAT(' Cannot find stretching factor', /, ' stretch = ', e15.6, &
           ' new_stretch = ', e15.6, ' iter = ', i3)
9100 FORMAT(' Stretching factor = ', e15.6, /, ' Match point       = ', &
           e15.6, /, ' First z           = ', e15.6, /, ' Number of iters   = ', i4)
5300 FORMAT(' nnz = ', i4)
5600 FORMAT(' 5600 in vgrid ', /, ' iz ', 5x, ' zw ', /, (i3, e15.6))

END SUBROUTINE
