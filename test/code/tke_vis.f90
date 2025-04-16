subroutine tke_vis(istep)
!GET VISCOSITY USING DEARDORFF RKE MODEL WITH STABILITY CORRECTION.
!FIXES FOR SURFACE LAYER
!GET RHS OF EQUATION
!only called in solve/comp1

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  REAL fnt1(nnx,iys:iye), fnt2(nnx,iys:iye,izs:ize)
  REAL fnt3(nnx,iys:iye)
  REAL ex(nnx,iys:iye), ey(nnx,iys:iye,izs:ize)
  REAL u_avg(nnx,iys:iye), v_avg(nnx,iys:iye), dissp(nnx,iys:iye)
  REAL alk(nnx,iys:iye,izs-1:ize+1)

  !GET LENGTH SCALES AND EDDY VISCOSITY
  IF(i_dear == 0) THEN !(i_dear=0 means deardorff, if /=0, then Schumann)
    CALL dear_vis(alk)
  ELSE
    CALL schu_vis(alk)
  ENDIF

  !IF SPECIAL 2 PART SURFACE LAYER MODEL IS ON
  !GET 'MEAN' VISCOSITY
  DO iz=izs-1,ize
    izm1         = iz - 1
    izp1         = iz + 1
    vis_mean(iz) = 0.0
    IF(ivis == 1 .AND. iz <= nmatch) THEN !ivis0 == 1 (new eddy viscosity model), then ivis=1
      IF(iz <= 1) THEN !if bottom boundary
        vis_mean(iz) = xksurf !xksurf defined in boundary/surfvis
      ELSE !it is not the bottom boundary 
        stravg = SQRT((u_mn(izp1)-u_mn(iz))**2 + (v_mn(izp1)-v_mn(iz))**2)* &
              ABS(dzu_i(izp1))
        vis_mean(iz) = xksurf*viscon*stravg
      ENDIF
    ENDIF
  ENDDO

  !UPDATE RHS OF SGS E FROM X AND Z PIECES
  !CUBE OF SIZE (NNX, IYZ, IYE, IZS:IZE)
  DO iz=izs,ize
    izm1   = iz - 1
    izp1   = iz + 1
    weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1  = 1.0 - weit
    dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
    dzw3_i = 2.0*dzw2_i
    dslk   = dsl_z(iz)

    DO iy=iys,iye
      DO ix=1,nnx
        ex(ix,iy) = e(ix,iy,iz)
      ENDDO
    ENDDO

    CALL xderivp(ex(1,iys),trigx(1,1),xk(1),nnx,iys,iye)

    !INCLUDE STOKES CONTRIBUTION IN ADVECTION AND HORIZONTAL X-DIFFUSION
    DO iy=iys,iye
      DO ix=1,nnx
        u_avg(ix,iy)  = (stokes(iz)*dir_x + u(ix,iy,iz))*weit1 +            &
              (stokes(izp1)*dir_x + u(ix,iy,izp1))*weit
        fnt1(ix,iy)   = e(ix,iy,iz)*u_avg(ix,iy) - 4.0*vis_m(ix,iy,iz)*ex(ix,iy)
      ENDDO
    ENDDO

    CALL xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)

    DO iy=iys,iye
      DO ix=1,nnx
        r5(ix,iy,iz) = -fnt1(ix,iy) - (w(ix,iy,izp1)*e(ix,iy,izp1) -        &
              w(ix,iy,izm1)*e(ix,iy,izm1))*dzw2_i
        r5(ix,iy,iz)=0.25*((r5(ix,iy,iz) - u_avg(ix,iy)*ex(ix,iy))*2.0      &
              - w(ix,iy,iz)*(e(ix,iy,izp1)-e(ix,iy,izm1))*dzw3_i)
      ENDDO
    ENDDO

    !9/1989 ADD ihflt=1 OPTION -- MEAN SHEAR DOES NOT GENERATE SGS TKE
    uxymm=0.
    uxymp=0.
    vxymm=0.
    vxymp=0.

    IF(ivis == 1 .AND. iz <= nmatch) THEN
      uxymm = u_mn(iz)
      uxymp = u_mn(izp1)
      vxymm = v_mn(iz)
      vxymp = v_mn(izp1)
    ENDIF

    DO iy=iys,iye
      DO ix=1,nnx
        !DISSIPATION
        dissp(ix,iy) =  (0.19+0.74*alk(ix,iy,iz)/dslk)*e(ix,iy,iz)*         &
              SQRT(e(ix,iy,iz))/alk(ix,iy,iz)
        r5(ix,iy,iz)=r5(ix,iy,iz) - dissp(ix,iy)

        !VERTICAL DIFFUSION
        fnt3(ix,iy) = ((vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*(e(ix,iy,izp1)-  &
              e(ix,iy,iz))*dzw_i(izp1) - (vis_m(ix,iy,iz)+                  &
              vis_m(ix,iy,izm1))*(e(ix,iy,iz  )-e(ix,iy,izm1))*dzw_i(iz))*  &
              dzu_i(izp1)
        r5(ix,iy,iz) = r5(ix,iy,iz) + fnt3(ix,iy)

        !SHEAR PRODUCTION
        s11 = weit1*ux(ix,iy,iz)**2 + weit*ux(ix,iy,izp1)**2
        s22 = weit1*vy(ix,iy,iz)**2 + weit*vy(ix,iy,izp1)**2
        wz  = (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
        wzp = (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
        s33 = weit*wzp**2 + weit1*wz**2
        s12 = weit1*(uy(ix,iy,iz) + vx(ix,iy,iz))**2 + weit*(uy(ix,iy,izp1) &
              + vx(ix,iy,izp1))**2
        uzmn=(u(ix,iy,izp1)-uxymp-u(ix,iy,iz)+uxymm)*dzu_i(izp1)
        vzmn=(v(ix,iy,izp1)-vxymp-v(ix,iy,iz)+vxymm)*dzu_i(izp1)
        s13 = (uzmn + wx(ix,iy,iz))**2
        s23 = (vzmn + wy(ix,iy,iz))**2

        fnt1(ix,iy) = vis_m(ix,iy,iz)*(2.0*(s11 + s22 + s33) + s13 + s23 + s12)
        r5(ix,iy,iz) = r5(ix,iy,iz) + fnt1(ix,iy)

        !GENERAL STOKES PRODUCTION
        dstdz   = (stokes(izp1) - stokes(iz))*dzu_i(izp1)
        st_prod = dstdz*dir_x*vis_m(ix,iy,iz)*(wx(ix,iy,iz) + uzmn) +       &
              dstdz*dir_y*vis_m(ix,iy,iz)*(wy(ix,iy,iz) + vzmn)
        fnt1(ix,iy) = fnt1(ix,iy) + st_prod
        r5(ix,iy,iz) = r5(ix,iy,iz)+fnt1(ix,iy)

        !BUOYANCY, GET TAU_W*THETA
        buoy_sgs = -vis_sv(ix,iy,iz)*(t(ix,iy,1,izp1)-t(ix,iy,1,iz))*dzu_i(izp1)
        r5(ix,iy,iz) = r5(ix,iy,iz) + batag*buoy_sgs
      ENDDO
    ENDDO

    !COMPUTE SHEAR, BUOYANCY, DIFFUSION
    !TERMS IN SGS E EQUATION FOR PRINTOUT
    !*** TRIZ IS ONLY IN VERTICAL DIFFUSION ***
    IF(istep == 1) THEN
      shrz(iz)   = 0.0
      triz(iz)   = 0.0
      t_diss(iz) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          shrz(iz)   = shrz(iz) + fnt1(ix,iy)
          t_diss(iz) = t_diss(iz) + dissp(ix,iy)
          triz(iz)   = triz(iz) + fnt3(ix,iy)
        ENDDO
      ENDDO
      shrz(iz)   = shrz(iz)*fnxy
      t_diss(iz) = t_diss(iz)*fnxy
      triz(iz)   = triz(iz)*fnxy
    ENDIF
  ENDDO

  !UPDATE TENDENCY OF SGS E FROM Y CONTRIBUTIONS
  !PENCIL SIZE (NNX, IYS:IYE, IZS:IZE)
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        ey(ix,iy,iz) = e(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  CALL yd_mpi(ey(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,iys, &
        iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

  !SKEW SYMMETRIC ADVEXTION [VDE/DY + D/DY(VE)]/2
  !PLUS SGS DIFFUSION CONTRIBUTION 
  DO iz=izs,ize
    izm1   = iz - 1
    izp1   = iz + 1
    weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1  = 1.0 - weit
    DO iy=iys,iye
      DO ix=1,nnx
        v_avg(ix,iy)   = (stokes(iz)*dir_y + v(ix,iy,iz))*weit1 +           &
              (stokes(izp1)*dir_y + v(ix,iy,izp1))*weit
        fnt2(ix,iy,iz) = e(ix,iy,iz)*v_avg(ix,iy) - 4.0*vis_m(ix,iy,iz)*    &
              ey(ix,iy,iz)
        r5(ix,iy,iz)   = r5(ix,iy,iz) - 0.5*(v_avg(ix,iy)*ey(ix,iy,iz))
      ENDDO
    ENDDO
  ENDDO

  CALL yd_mpi(fnt2(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,   &
        iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        r5(ix,iy,iz) = r5(ix,iy,iz) - 0.5*fnt2(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
