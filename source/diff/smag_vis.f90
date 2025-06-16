SUBROUTINE smag_vis
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  implicit none

  REAL :: sij2
  REAL :: alk(nnx, iys:iye, izs-1:ize+1)
  real :: dslk, s11, s22, s33, s12, s13, s23, wz, wzp, weit, weit1, uzmn, vzmn
  integer :: i, j, iz

  DO iz=izs-1,ize+1
    vis_mean(iz) = 0.0
    dslk  = dsl_z(iz)
    IF(iz .gt. 0) dslk  = AMIN1(dsl_z(iz),vk*ABS(z(iz))/csmag)

    alk(:, :, iz) = dslk

    weit   = dzw(iz)/(dzw(iz) + dzw(iz+1))
    weit1  = 1.0 - weit
    !horizontal average of fluctuations
    !ufluct_vfluct_mn = 0.0
    !vfluct_Tfluct_mn = 0.0
    !DO j=iys,iye
    !  DO i=1,nnx
    !    ufluct_vfluct_mn = ufluct_vfluct_mn + (u(i, j, iz+1) - u_mn(iz)) * (v(i, j, iz+1) - v_mn(iz))
    !    vfluct_Tfluct_mn = vfluct_Tfluct_mn + (v(i, j, iz+1) - v_mn(iz)) * (t(i, j, iz+1) - t_mn(iz))
    !  END DO
    !ENDDO
    !ufluct_vfluct_mn = ufluct_vfluct_mn*fnxy
    !vfluct_Tfluct_mn = vfluct_Tfluct_mn*fnxy
    ! GET STRAINS and compute viscosity
    DO j=iys,iye
      DO i=1,nnx
        s11 = weit1*ux(i, j, iz)**2 + weit*ux(i, j, iz+1)**2
        s22 = weit1*vy(i, j, iz)**2 + weit*vy(i, j, iz+1)**2
        wz = (w(i, j, iz)-w(i, j, iz-1))*dzw_i(iz)
        wzp = (w(i, j, iz+1)-w(i, j, iz))*dzw_i(iz+1)
        s33 = weit*wzp**2 + weit1*wz**2
        s12 = weit1*(uy(i, j, iz) + vx(i, j, iz))**2 + weit*(uy(i, j, iz+1) &
              + vx(i, j, iz+1))**2
        uzmn=(u(i, j, iz+1)-u(i, j, iz))*dzu_i(iz+1)
        vzmn=(v(i, j, iz+1)-v(i, j, iz))*dzu_i(iz+1)
        s13 = (uzmn + wx(i, j, iz))**2
        s23 = (vzmn + wy(i, j, iz))**2

        sij2 = (s11 + s22 + s33) + 0.5 * (s13 + s23 + s12)
        !eddy viscosity
        vis_m(i, j, iz)  = (csmag*dslk)**2 * SQRT(2 * sij2)

        !eddy diffusivity
        !dudy_mean = (u_mn(iz+1)-u_mn(iz)) / dy
        !dTdy_mean = (t_mn(iz+1)-t_mn(iz)) / dy
        vis_s(i, j, iz)  = vis_m(i, j, iz) !* ufluct_vfluct_mn * dTdy_mean / (vfluct_Tfluct_mn * dudy_mean)
        vis_sv(i, j, iz) = vis_s(i, j, iz)
      ENDDO
    ENDDO

    ! SPECIAL CASE FOR IZ = 1
    IF(iz==1 .AND. ibcl == 0) THEN
      DO j=iys,iye
        DO i=1,nnx
          vis_m(i, j, iz-1)  = vis_m(i, j, iz)
          vis_s(i, j, iz-1)  = vis_s(i, j, iz)
          vis_sv(i, j, iz-1) = vis_sv(i, j, iz)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE
