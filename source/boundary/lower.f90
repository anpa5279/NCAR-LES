SUBROUTINE lower(it)
! SETUP LOWER BC FOR PLANE AT IZ = 1
! CAN USE BUSINGER OR LARGE FORMULAS WITH WIND
! INDEX F(.,.,2) INDICATES LOWER.
! THREADED VERSION

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  REAL :: sfc_flx(2+nscl)

  iz   = 1
  izm1 = iz - 1
  dz_i = dzu_i(1)

  DO iy=iys,iye
    DO ix=1,nnx
      ebc(ix,iy,2)  = AMAX1(e(ix,iy,iz),sml_eg)
      wbc(ix,iy,2)  = 0.0
      pbc(ix,iy,2)  = 0.0
      pbc2(ix,iy,2) = 0.0
    ENDDO
  ENDDO

  CALL sufto(it)

  DO iy=iys,iye
    DO ix=1,nnx
      tau13m(ix,iy) = -au13m
      tau23m(ix,iy) = -au23m
    ENDDO
  ENDDO

  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        taut3m(ix,iy,iscl) = aut3m(iscl)
      ENDDO
    ENDDO
  ENDDO

  ! PARTIAL SUMS OF SURFACE FLUXES AND MEAN SCALAR
  sfc_flx(1) = 0.0
  sfc_flx(2) = 0.0
  DO iy=iys,iye
    DO ix=1,nnx
      sfc_flx(1) = sfc_flx(1) + tau13m(ix,iy)
      sfc_flx(2) = sfc_flx(2) + tau23m(ix,iy)
    ENDDO
  ENDDO

  DO iscl=1,nscl
    sfc_flx(2+iscl) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        sfc_flx(2+iscl) = sfc_flx(2+iscl) + taut3m(ix,iy,iscl)
      ENDDO
    ENDDO
  ENDDO

  CALL mpi_sum_xy(sfc_flx,myid,iss,ise,(2+nscl))

  uwsfc = sfc_flx(1)*fnxy
  vwsfc = sfc_flx(2)*fnxy

  DO iscl=1,nscl
    wtsfc(iscl) = sfc_flx(2+iscl)*fnxy
  ENDDO

  WRITE(nprt,2345) uwsfc, vwsfc, wtsfc(nscl), tsfcc(nscl)

  DO iy=iys,iye
    DO ix=1,nnx
      dudz     = 2.*(u(ix,iy,iz) + ugal)*dz_i
      dvdz     = 2.*v(ix,iy,iz)*dz_i
      ubc(ix,iy,2) = u(ix,iy,iz) !- dudz*dzu(iz)
      vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
    ENDDO
  ENDDO

  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
        tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
      ENDDO
    ENDDO
  ENDDO

  ! INITIALIZE U,V,W,T AND DERIVATIVES AT IZM1
  DO iy=iys,iye
    DO ix=1,nnx
      u(ix,iy,izm1)  = ubc(ix,iy,2)
      v(ix,iy,izm1)  = vbc(ix,iy,2)
      w(ix,iy,izm1)  = wbc(ix,iy,2)
      r3(ix,iy,izm1) = 0.0
      e(ix,iy,izm1)  = ebc(ix,iy,2)
      ux(ix,iy,izm1) = 0.0
      uy(ix,iy,izm1) = 0.0
      vx(ix,iy,izm1) = 0.0
      vy(ix,iy,izm1) = 0.0
      wx(ix,iy,izm1) = wbc(ix,iy,2)
      wy(ix,iy,izm1) = wbc(ix,iy,2)
    ENDDO
  ENDDO

  ! NO NEED TO CALL DERIVATIVES HERE SINCE WBC = 0
  ! CHANGE FOR MORE GENERAL LOWER BC
  DO iscl=1,nscl
    DO iy=iys,iye
      DO ix=1,nnx
        t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
      ENDDO
    ENDDO
  ENDDO

  RETURN

! FORMAT
2345  FORMAT(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,' wtsfc = ', &
            e15.6,' tsfcc = ',e15.6)

END SUBROUTINE
