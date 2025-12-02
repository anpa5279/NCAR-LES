SUBROUTINE solve_trid(pt, ptop)
! TRIDIAGONAL SOLVER
! ODD ORDER FOR PTOP, PTOP2, BECAYSE OF 2D FFT

  USE pars
  USE con_data
  USE con_stats

  REAL :: ptop(nny,jxs:jxe,1:2)
  REAL :: pt(0:nnz+1,jxs:jxe,jys:jye)
  REAL :: aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
  REAL :: fac_u(nnz), fac_l(nnz), fac_a(nnz)

  DO iz=1,nnz
    fac_u(iz) = 1.0/(dzw(iz)*dzu(iz+1))
    fac_l(iz) = 1.0/(dzw(iz)*dzu(iz))
    fac_a(iz) = fac_l(iz) + fac_u(iz)
  ENDDO

  DO kp=jys,jye
    DO lp=jxs,jxe
      DO iz=2,nnz-1
        bb(iz,lp)  = fac_l(iz)
        aa(iz,lp)  = fac_u(iz)
        dd(iz,lp)  = -xks(lp,kp) - fac_a(iz)
        rh(iz,lp)  = pt(iz,lp,kp)
      ENDDO
    ENDDO

    ! LOWER BOUNDARY, FILL EXTERIOR PRESSURE (NOT USED)
    DO lp=jxs,jxe
      bb(1,lp)  = 1.0
      aa(1,lp)  = fac_u(1)
      dd(1,lp)  = -xks(lp,kp) - fac_u(1)
      rh(1,lp)  = pt(1,lp,kp)
      pt(0,lp,kp) = 0.0
    ENDDO

    ! UPPER BOUNDARY, FILL EXTERIOR PRESSURE (NOT USED)
    IF(ibcu == 1) THEN
      DO lp=jxs,jxe
        bb(nnz,lp) = 0.0
        aa(nnz,lp) = 0.0
        dd(nnz,lp) = 1.0
        rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
        pt(nnz+1,lp,kp) = 0.0
      ENDDO
    ELSE
      DO lp=jxs,jxe
        bb(nnz,lp) = fac_l(nnz)
        aa(nnz,lp) = 1.0
        dd(nnz,lp) = -xks(lp,kp) - fac_l(nnz)
        rh(nnz,lp) = pt(nnz,lp,kp)
        pt(nnz+1,lp,kp) = 0.0
      ENDDO
    ENDIF

    ! SPECIAL SITUATION FOR ZEROTH MODE MAKES MEAN PRESSURE = 0
    IF(kp == 1 .AND. jxs == 1) THEN
      DO iz=1,nnz
        dd(iz,1) = 1.0
        rh(iz,1) = 0.0
        aa(iz,1) = 0.0
        bb(iz,1) = 0.0
        dd(iz,2) = 1.0
        rh(iz,2) = 0.0
        aa(iz,2) = 0.0
        bb(iz,2) = 0.0
      ENDDO
    ENDIF

    ! SOLVE SYSTEM
    CALL tridv(bb,dd,aa,rh,nnz,jxs,jxe)

    DO lp=jxs,jxe
      DO iz=1,nnz
        pt(iz,lp,kp) = rh(iz,lp)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
