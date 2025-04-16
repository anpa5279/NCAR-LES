SUBROUTINE rhs_uvw(istep)
!GET RHS OF (U,V,W) EQUATIONS FOR PENCIL SIZE (NNX,IYS:IYE,IZS:IZE)
!called in solve/comp1 only
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  REAL :: fntd(nnx,iys:iye,izs:ize)
  REAL :: fnt1(nnx,iys:iye), fnt2(nnx,iys:iye)
  REAL :: fnt3(nnx,iys:iye), fnt4(nnx,iys:iye)
  REAL :: tau13_u(nnx,iys:iye), tau23_u(nnx,iys:iye)
  REAL :: tau13_l(nnx,iys:iye), tau23_l(nnx,iys:iye)
  REAL :: r3_sum(1:nnz)

  DO iz=izs,ize
    izm1 = iz - 1
    izp1 = iz + 1
    weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1 = 1.0 - weit

    !DYNAMICS
    DO iy=iys,iye
      DO ix=1,nnx
        uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
        vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
        uz  = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
        vz  = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)

        u_avg = u(ix,iy,iz)*weit1 + u(ix,iy,izp1)*weit
        v_avg = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit

        !ADVECTION
        u_adv =  v(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz)) - 0.5*(w(ix,iy,iz) &
              *(uz - wx(ix,iy,iz))+w(ix,iy,izm1)*(uzm - wx(ix,iy,izm1)))
        v_adv = -u(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz)) + 0.5*(w(ix,iy,iz) &
              *(wy(ix,iy,iz) - vz)+w(ix,iy,izm1)*(wy(ix,iy,izm1) - vzm))
        w_adv = u_avg*(uz - wx(ix,iy,iz)) - v_avg*(wy(ix,iy,iz) - vz)

        !CORIOLIS, VERTICAL, AND HORIZONTAL COMPONENTS
        u_cor =  fcor*(v(ix,iy,iz) + stokes(iz)*dir_y) - fcor_h*w(ix,iy,iz)
        v_cor = -fcor*(u(ix,iy,iz) + stokes(iz)*dir_x)
        w_cor =  fcor_h*u(ix,iy,iz)

        !BUOYANCY (WITH HYDROSTATIC PART)
        w_buy = batag*(t(ix,iy,1,iz)*weit1 + t(ix,iy,1,izp1)*weit)

        !GEOSTROPHIC WIND
        u_geo = -fcor*vg(iz)
        v_geo =  fcor*(ug(iz)-ugal)

        !TOTALS
        r1(ix,iy,iz) = u_adv + u_cor + u_geo
        r2(ix,iy,iz) = v_adv + v_cor + v_geo
        r3(ix,iy,iz) = w_adv + w_cor + w_buy
      ENDDO
    ENDDO

    !STOKES TERM
    stokavg = stokes(iz)*weit1 + stokes(izp1)*weit
    DO iy=iys,iye
      DO ix=1,nnx
        r1(ix,iy,iz) = r1(ix,iy,iz) + stokes(iz)*dir_y*(vx(ix,iy,iz) -      &
              uy(ix,iy,iz))
        r2(ix,iy,iz) = r2(ix,iy,iz) + stokes(iz)*dir_x*(uy(ix,iy,iz) -      &
              vx(ix,iy,iz))
        uz = (u(ix,iy,izp1) - u(ix,iy,iz))*dzu_i(izp1)
        vz = (v(ix,iy,izp1) - v(ix,iy,iz))*dzu_i(izp1)
        r3(ix,iy,iz) = r3(ix,iy,iz) + stokavg*dir_x*(uz - wx(ix,iy,iz)) -   &
              stokavg*dir_y*(wy(ix,iy,iz) - vz)
      ENDDO
    ENDDO

    !GET TAU_13,_23 AT IZ-1
    IF (iz/=1 .OR. ibcl/=0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
          vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
          tau13_l(ix,iy) = -vis_m(ix,iy,izm1)*(uzm + wx(ix,iy,izm1)) -      &
                vis_mean(izm1)*(u_mn(iz)-u_mn(izm1))*dzu_i(iz)
          tau23_l(ix,iy) = -vis_m(ix,iy,izm1)*(vzm + wy(ix,iy,izm1)) -      &
                vis_mean(izm1)*(v_mn(iz)-v_mn(izm1))*dzu_i(iz)
        ENDDO
      ENDDO
    ELSE
      DO iy=iys,iye
        DO ix=1,nnx
          tau13_l(ix,iy) = tau13m(ix,iy)
          tau23_l(ix,iy) = tau23m(ix,iy)
        ENDDO
      ENDDO
    ENDIF

    !X AND Z HORIZONTAL SGS FLUXES FOR U,V,W
    !TAU_11,_12,_13,_23 AT IZ
    DO iy=iys,iye
      DO ix=1,nnx
        fnt1(ix,iy) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*ux(ix,iy,iz)
        fnt2(ix,iy) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*(uy(ix,iy,iz) &
              +vx(ix,iy,iz))
        uz = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
        vz = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
        tau13_u(ix,iy) = -vis_m(ix,iy,iz)*(uz+wx(ix,iy,iz)) - vis_mean(iz)* &
              (u_mn(izp1)-u_mn(iz))*dzu_i(izp1)
        tau23_u(ix,iy) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -             &
              vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
        fnt3(ix,iy) = tau13_u(ix,iy)
      ENDDO
    ENDDO

    CALL xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
    CALL xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
    CALL xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)

    DO iy=iys,iye
      DO ix=1,nnx
        r1(ix,iy,iz) = r1(ix,iy,iz) - fnt1(ix,iy)-(tau13_u(ix,iy)-          &
              tau13_l(ix,iy))*dzw_i(iz)
        r2(ix,iy,iz) = r2(ix,iy,iz) - fnt2(ix,iy)-(tau23_u(ix,iy)-          &
              tau23_l(ix,iy))*dzw_i(iz)
        fnt4(ix,iy) = -(vis_m(ix,iy,izm1)+vis_m(ix,iy,iz))*(w(ix,iy,iz)-    &
              w(ix,iy,izm1))*dzw_i(iz)
        fnt2(ix,iy) = -(vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*(w(ix,iy,izp1)-  &
              w(ix,iy,iz))*dzw_i(izp1)
        r3(ix,iy,iz) = r3(ix,iy,iz) - fnt3(ix,iy) - (fnt2(ix,iy)-           &
              fnt4(ix,iy))*dzu_i(izp1)
      ENDDO
    ENDDO

    !SAVE SGS FLUXES FOR PRINTOUT
    IF(istep == 1) THEN
      uwsb(iz)   = 0.0
      vwsb(iz)   = 0.0
      wwsb(iz)   = 0.0
      tr_tau(iz) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          uwsb(iz) = uwsb(iz) + tau13_u(ix,iy)
          vwsb(iz) = vwsb(iz) + tau23_u(ix,iy)
          wwsb(iz) = wwsb(iz) + fnt4(ix,iy)
          ufluc    = (u(ix,iy,izp1) - uxym(izp1))*weit + (u(ix,iy,iz) -     &
                uxym(iz))*weit1
          vfluc    = (v(ix,iy,izp1) - vxym(izp1))*weit + (v(ix,iy,iz) -     &
                vxym(iz))*weit1
          tr_tau(iz) = tr_tau(iz) + tau13_u(ix,iy)*ufluc + tau23_u(ix,iy)*vfluc
        ENDDO
      ENDDO
      uwsb(iz)   = uwsb(iz)*fnxy
      vwsb(iz)   = vwsb(iz)*fnxy
      wwsb(iz)   = wwsb(iz)*fnxy
      tr_tau(iz) = tr_tau(iz)*fnxy
    ENDIF
  ENDDO

  !SGS FLUXES TAU_12,_22,_23 THAT DEPEND ON Y-DERIVS
  DO iz=izs,ize
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx
        fntd(ix,iy,iz) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*           &
              (uy(ix,iy,iz)+vx(ix,iy,iz))
      ENDDO
    ENDDO
  ENDDO

  CALL yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,   &
        iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

  DO iz=izs,ize
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx
        r1(ix,iy,iz)   = r1(ix,iy,iz) - fntd(ix,iy,iz)
        fntd(ix,iy,iz) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*vy(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  CALL yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,   &
        iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

  DO iz=izs,ize
    izp1 = iz + 1
    DO iy=iys,iye
      DO ix=1,nnx
        r2(ix,iy,iz)   = r2(ix,iy,iz) - fntd(ix,iy,iz)
        vz             = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
        fntd(ix,iy,iz) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -             &
              vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
      ENDDO
    ENDDO
  ENDDO

  CALL yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,   &
        iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

  DO iz=1,nnz
    r3_sum(iz) = 0.0
  ENDDO

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        r3(ix,iy,iz) = r3(ix,iy,iz) - fntd(ix,iy,iz)
        r3_sum(iz)   = r3_sum(iz) + r3(ix,iy,iz)
      ENDDO
    ENDDO
    r3_sum(iz) = r3_sum(iz)*fnxy
  ENDDO

  CALL mpi_sum_z(r3_sum,i_root,myid,nnz,1)

  !MAKE SURE <R3> = 0 AND SET R3 = 0 AT TOP
  DO iz=izs,ize
    IF(iz == nnz) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          r3(ix,iy,iz) = 0.0
        ENDDO
      ENDDO
    ELSE
      DO iy=iys,iye
        DO ix=1,nnx
          r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  RETURN
END
