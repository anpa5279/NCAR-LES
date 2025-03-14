SUBROUTINE solve_trid(pt, ptop)
!TRIDIAGONAL SOLVER
!ODD ORDER FOR PTOP, PTOP2, BECAYSE OF 2D FFT
!only called in solve/pressure with variables (pt=pt, ptop=ptopfft)

  USE pars
  USE con_data
  USE con_stats

  REAL :: ptop(nny,jxs:jxe,1:2) !only used in line 47 when upper boundary is set by radiation (ibcu==1)
  REAL :: pt(0:nnz+1,jxs:jxe,jys:jye) !all z interfaces, starting and ending process lengths for x and y
          !aa= above diagonal, bb= below diagonal, dd=diagonal, rh= RHS
  REAL :: aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),dd(nnz,jxs:jxe),rh(nnz,jxs:jxe) !number of layers in z, starting and ending process lengths for x and y
  REAL :: fac_u(nnz), fac_l(nnz), fac_a(nnz) !lnumber of layers in z

  DO iz=1,nnz
    fac_u(iz) = 1.0/(dzw(iz)*dzu(iz+1))
    fac_l(iz) = 1.0/(dzw(iz)*dzu(iz))
    fac_a(iz) = fac_l(iz) + fac_u(iz)
  ENDDO

  DO kp=jys,jye !starting process and ending process for y. defined in setup/set_range
    DO lp=jxs,jxe !starting process and ending process for x. defined in setup/set_range
      DO iz=2,nnz-1 !looking at inner meshing only
        !saving to z-y matrix
        bb(iz,lp)  = fac_l(iz)
        aa(iz,lp)  = fac_u(iz)
        dd(iz,lp)  = -xks(lp,kp) - fac_a(iz)
        rh(iz,lp)  = pt(iz,lp,kp)
      ENDDO
    ENDDO

    !LOWER BOUNDARY, FILL EXTERIOR PRESSURE (NOT USED)
    DO lp=jxs,jxe !bottom boundary mesh 
      bb(1,lp)  = 1.0
      aa(1,lp)  = fac_u(1)
      dd(1,lp)  = -xks(lp,kp) - fac_u(1)
      rh(1,lp)  = pt(1,lp,kp)
      pt(0,lp,kp) = 0.0
    ENDDO

    !UPPER BOUNDARY, FILL EXTERIOR PRESSURE (NOT USED)
    IF(ibcu == 1) THEN !ibcu    = 1 ; upper boundary condition set by radiation bc
      DO lp=jxs,jxe !applying conditions to the top boundary
        bb(nnz,lp) = 0.0
        aa(nnz,lp) = 0.0
        dd(nnz,lp) = 1.0
        rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
        pt(nnz+1,lp,kp) = 0.0 !top interface
      ENDDO
    ELSE !ibcu (=0 or -1) is a fixed value of 0 or is defined by a coarse mesh
      DO lp=jxs,jxe !applying conditions to the top boundary
        bb(nnz,lp) = fac_l(nnz)
        aa(nnz,lp) = 1.0
        dd(nnz,lp) = -xks(lp,kp) - fac_l(nnz)
        rh(nnz,lp) = pt(nnz,lp,kp)
        pt(nnz+1,lp,kp) = 0.0
      ENDDO
    ENDIF

    !SPECIAL SITUATION FOR ZEROTH MODE MAKES MEAN PRESSURE = 0
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

    !SOLVE SYSTEM
    CALL tridv(bb,dd,aa,rh,nnz,jxs,jxe) !TRIDIAGONAL MATRIX SOLVER WITH MULTIPLE VECTORS. output rh

    DO lp=jxs,jxe
      DO iz=1,nnz
        pt(iz,lp,kp) = rh(iz,lp)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
