SUBROUTINE smag_vis(istep)
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: d_grid(izs-1:ize+1)
  REAL :: sij2(nnx,iys:iye,izs-1:ize+1)

  DO iz=izs,MIN(ize,nmatch)
    izp1 = iz + 1
    izm1 = iz - 1
    d_grid(iz) = (ABS(dx*dy*dzw(iz)))**(1./3.)
    ! GET FLUCTUATING STRAINS
    DO j=iys,iye
      DO i=1,nnx
        s11 = (ux(i,j,iz))**2
        s22 = (vy(i,j,iz))**2
        wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
        s33 =  (wzp)**2
        s12 =  (0.5 * (uy(i,j,iz) + vx(i,j,iz)))**2
        s13 =  (0.5 * ((((u(i,j,izp1) - u(i,j,iz)) * dzu_i(izp1) + wx(i,j,iz)))))**2
        s23 =  (0.5 * ((v(i,j,izp1) - v(i,j,iz)) * dzu_i(izp1) + wy(i,j,iz)))**2
        sij2(i,j,iz) = s11 + s22 + s33 + 2 * s12 + 2 * s13 + 2 * s23
        vis_m(i,j,iz)  = (csmag)**2 * (d_grid(iz))**2 * SQRT(2 * sij2(i, j, iz))
        vis_s(i,j,iz)  = vis_m(i,j,iz)
        vis_sv(i,j,iz) = vis_s(i,j,iz)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE