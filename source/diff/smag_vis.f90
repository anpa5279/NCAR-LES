SUBROUTINE smag_vis(istep)
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: d_grid(izs-1:ize+1)
  REAL :: sij2(nnx,iys:iye,izs-1:ize+1)

  DO iz=izs-1,ize+1
    izm1   = iz - 1
    izp1 = iz + 1
    weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
    weit1  = 1.0 - weit
    d_grid(iz) = (ABS(dx*dy*dzw(iz)))**(1./3.)
    DO j=iys,iye
      DO i=1,nnx
        !strain rate tensor terms squared. strain rate tens: sij= 1/2*(u_<i,j> + u_<j,i>) 
        s11 = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
        s22 = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
        wz  = (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
        wzp = (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
        s33 = (weit*wzp + weit1*wz)**2
        s12 = weit1*(uy(ix,iy,iz) + vx(ix,iy,iz))**2 + weit*(uy(ix,iy,izp1) &
              + vx(ix,iy,izp1))**2
        uzmn=(u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
        vzmn=(v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
        s13 = (uzmn + wx(ix,iy,iz))**2
        s23 = (vzmn + wy(ix,iy,iz))**2
        !s12 =  (0.5 * (weit1*(uy(i,j,iz) + vx(i,j,iz)) + weit*(uy(i,j,izp1) + vx(i,j,izp1))))**2
        !s13 =  (0.5 * ((((u(i,j,izp1) - u(i,j,iz)) * dzu_i(izp1) + wx(i,j,iz)))))**2
        !s23 =  (0.5 * ((v(i,j,izp1) - v(i,j,iz)) * dzu_i(izp1) + wy(i,j,iz)))**2
        sij2(i,j,iz) = 2.0*(s11 + s22 + s33) + s13 + s23 + s12!s11 + s22 + s33 + 2.0 * s12 + 2.0 * s13 + 2.0 * s23
        vis_m(i,j,iz)  = (csmag)**2 * (d_grid(iz))**2 * SQRT(2.0 * sij2(i, j, iz))
        vis_s(i,j,iz)  = vis_m(i,j,iz)
        vis_sv(i,j,iz) = vis_s(i,j,iz)
        r5(ix,iy,iz) = 0.0
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE