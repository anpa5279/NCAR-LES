SUBROUTINE smag_vis
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: sij2(nnx,iys:iye,izs-1:ize+1)
  REAL alk(nnx,iys:iye,izs-1:ize+1)
  REAL :: dslk, weit, weit1
  REAL :: s11, s22, s33, s12, s23, s13

  DO iz=izs-1,ize+1
    izp1 = iz + 1
    dslk  = dsl_z(iz)
    IF(iz .gt. 0) dslk  = AMIN1(dsl_z(iz),vk*ABS(z(iz))/csmag)
    ! NO STABILITY CORRECTED LENGTH SCALES
    DO j=iys,iye
      DO i=1,nnx
        alk(i,j,iz) = dslk
      END DO
    END DO

    weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1  = 1.0 - weit
    !horizontal average of fluctuations
    !ufluct_vfluct_mn = 0.0
    !vfluct_Tfluct_mn = 0.0
    !DO j=iys,iye
    !  DO i=1,nnx
    !    ufluct_vfluct_mn = ufluct_vfluct_mn + (u(i,j,izp1) - u_mn(iz)) * (v(i,j,izp1) - v_mn(iz))
    !    vfluct_Tfluct_mn = vfluct_Tfluct_mn + (v(i,j,izp1) - v_mn(iz)) * (t(i,j,izp1) - t_mn(iz))
    !  END DO
    !ENDDO
    !ufluct_vfluct_mn = ufluct_vfluct_mn*fnxy
    !vfluct_Tfluct_mn = vfluct_Tfluct_mn*fnxy
    ! GET STRAINS and compute viscosity
    DO j=iys,iye
      DO i=1,nnx
        s11 = weit1*ux(ix,iy,iz)**2 + weit*ux(ix,iy,izp1)**2
        s22 = weit1*vy(ix,iy,iz)**2 + weit*vy(ix,iy,izp1)**2
        wz  = (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
        wzp = (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
        s33 = weit*wzp**2 + weit1*wz**2
        s12 = weit1*(uy(ix,iy,iz) + vx(ix,iy,iz))**2 + weit*(uy(ix,iy,izp1) &
              + vx(ix,iy,izp1))**2
        uzmn=(u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
        vzmn=(v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
        s13 = (uzmn + wx(ix,iy,iz))**2
        s23 = (vzmn + wy(ix,iy,iz))**2

        sij2(ix,iy, iz) = (s11 + s22 + s33) + 0.5 * (s13 + s23 + s12)
        !eddy viscosity
        vis_m(i,j,iz)  = (csmag*dslk)**2 * SQRT(2 * sij2(i, j, iz))
        
        !eddy diffusivity
        !dudy_mean = (u_mn(izp1)-u_mn(iz)) / dy
        !dTdy_mean = (t_mn(izp1)-t_mn(iz)) / dy
        vis_s(i,j,iz)  = vis_m(i,j,iz) !* ufluct_vfluct_mn * dTdy_mean / (vfluct_Tfluct_mn * dudy_mean)
        vis_sv(i,j,iz) = vis_s(i,j,iz)
      ENDDO
    ENDDO
    print*, "iz", iz, "csmag", csmag! "dzw_i(izp1)", dzw_i(izp1), "dzu_i(izp1)", dzu_i(izp1)
    print*, "vis_m", vis_m(1,1,iz), "vis_s", vis_s(1,1,iz), "dslk", dslk, "sij2", sij2(1,1,iz)
    ! SPECIAL CASE FOR IZ = 1
    IF(iz==1 .AND. ibcl == 0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          vis_m(ix,iy,iz-1)  = vis_m(ix,iy,iz)
          vis_s(ix,iy,iz-1)  = vis_s(ix,iy,iz)
          vis_sv(ix,iy,iz-1) = vis_sv(ix,iy,iz)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE
