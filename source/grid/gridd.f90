SUBROUTINE gridd
! ALLOCATE SPACE AND PASS ARRAYS USING MODULES

    USE pars
    USE fields
    USE fftwk
    USE con_data
    USE con_stats

    ! ESTABLISH ASSOCIATION BETWEEN POINTERS AND DATA STRUCTURES
    CALL fill_cc
    CALL fill_cs

    ! DEBUG FOR ARRAYS
    big = -99.0e+300

    ! SETUP GRID
    nnx = nxg1
    nny = nyg1
    nnz = nzg1

    ! MAKE SURE PROBLEM CPU'S MATCH
    maxp = numprocs - 1
    ncpu_z = numprocs / ncpu_s
    IF (MOD(numprocs, ncpu_s) /= 0 .OR. ncpu_z > nnz) THEN
        IF (l_root) WRITE (6, 1000) numprocs, ncpu_s, mmz
        WRITE (nprt, 1000) numprocs, ncpu_s, nnz
        CALL mpi_finalize(ierr)
        STOP
    ELSE
        IF (l_root) WRITE (6, 1100) ncpu_s, ncpu_z, numprocs, maxp

        WRITE (nprt, 1100) ncpu_s, ncpu_z, numprocs, maxp

        ! ALLOCATE ARRAYS FOR (I,J,K) INDEXING ON EACH PROCESSOR (SEE SET_RANGE)
        ALLOCATE (ix_s(0:maxp), ix_e(0:maxp), jx_s(0:maxp), jx_e(0:maxp), &
                  kx_s(0:maxp), kx_e(0:maxp), mx_s(0:maxp), mx_e(0:maxp), &
                  iy_s(0:maxp), iy_e(0:maxp), jy_s(0:maxp), jy_e(0:maxp), &
                  is_s(0:maxp), is_e(0:maxp), iz_s(0:maxp), iz_e(0:maxp))

        ! SETUP ARRAY SIZES AND VARIABLE DIMENSION
        nxy = nnx * nny
        ncx = nnx / 2 + 1
        ncy = nny / 2 + 1
        nnxp1 = nnx + 1
        nnyp1 = nny + 1
        nnxp2 = nnx + 2
        nnyp2 = nny + 2
        nnzp1 = nnz + 1
        nnzm1 = nnz - 1
        ivis = ivis0
        fnxy = 1.0 / FLOAT(nnx * nny)

        WRITE (nprt, 7001) nnx, nny, nnz

        CALL set_range

        WRITE (nprt, 7002) nnx, nny, nnz

        num_y = iye + 1 - iys

        ! ALLOCATE SOLUTION ARRAYS
        ! ACCOUNT FOR NNXP2 FOR FIELDS BUT NOT ON RHS, POSSIBLE MONOTONE FOR SCALARS
        ALLOCATE (u(nnxp2, iys:iye, izs - 1:ize + 1), v(nnxp2, iys:iye, izs - 1:ize + 1), &
                  w(nnxp2, iys:iye, izs - 1:ize + 1), t(nnxp2, iys:iye, nscl, izs - 2:ize + 2), &
                  e(nnxp2, iys:iye, izs - 1:ize + 1), r1(nnx, iys:iye, izs - 1:ize + 1), &
                  r2(nnx, iys:iye, izs - 1:ize + 1), r3(nnx, iys:iye, izs - 1:ize + 1), &
                  r4(nnx, iys:iye, nscl, izs - 1:ize + 1), r5(nnx, iys:iye, izs - 1:ize + 1))

        ! ALLOCATE SPACE FOR BC ARRAYS ON TOP AND BOTTOM OF DOMAIN
        ALLOCATE (ubc(nnx, iys:iye, 2), vbc(nnx, iys:iye, 2), wbc(nnx, iys:iye, 2), &
                  tbc(nnx, iys:iye, nscl, 2), ebc(nnx, iys:iye, 2), pbc(nnx, iys:iye, 2), &
                  pbc2(nnx, iys:iye, 2))

        ! ALLOCATE SPACE FOR WIND AND SURFACE ARRAYS
        ALLOCATE (wind(nnx, iys:iye), tau13m(nnx, iys:iye), tau23m(nnx, iys:iye), &
                  taut3m(nnx, iys:iye, nscl), t_grnd(nnx, iys:iye, nscl))

        ! ALLOCATE SPACE FOR DERIVATIVE ARRAYS
        ALLOCATE (ux(nnx, iys:iye, izs - 1:ize + 1), uy(nnx, iys:iye, izs - 1:ize + 1), &
                  vx(nnx, iys:iye, izs - 1:ize + 1), vy(nnx, iys:iye, izs - 1:ize + 1), &
                  wx(nnx, iys:iye, izs - 1:ize + 1), wy(nnx, iys:iye, izs - 1:ize + 1))

        ! ALLOCATE SPACE FOR PRESSURE AND PRESSURE BCS
        ALLOCATE (p(nnxp2, iys:iye, izs - 1:ize + 1), ptop(nnxp2, iys:iye, 2))

        ! ALLOCATE SPACE FOR VISCOSITY AND DIFFUSIVITY
        ALLOCATE (vis_m(nnx, iys:iye, izs - 1:ize + 1), vis_s(nnx, iys:iye, izs - 1:ize + 1) &
                  , vis_sv(nnx, iys:iye, izs - 1:ize + 1))

        ! ALLOCATE SPACE FOR FFT TRIG FACTORS
        nq_trig = MAX(nnx, nny)
        ALLOCATE (trigx(2 * nq_trig + 15, 2), trigc(4 * nq_trig + 15))
    END IF

    RETURN

! FORMAT
1100 FORMAT(' Number of x-y slab cpus = ', i5, /, &
           ' Number of z-level cpus  = ', i5, /, &
           ' Total number of cpus    = ', i5, /, ' Max-p for index arrays  = ', i5)
7001 FORMAT(' 7001 gridd nnx = ', i4, ' nny = ', i4, ' nnz = ', i4)
7002 FORMAT(' 7002 gridd nnx = ', i4, ' nny = ', i4, ' nnz = ', i4)
1000 FORMAT(' Gridd Trouble number of processors and grid', &
           ' partitioning do not match.', /, ' Total num of cpus   = ', i5, &
           ' Num cpu on x-y slab = ', i5, /, ' Num of z-levels     = ', i5)

END SUBROUTINE