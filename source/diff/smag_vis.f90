SUBROUTINE smag_vis(alk)
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: sij2(nnx,iys:iye,izs-1:ize+1)
  REAL alk(nnx,iys:iye,izs-1:ize+1)
  DO iz=izs-1,ize+1
    izp1 = iz + 1
    izm1 = iz - 1
    dslk  = dsl_z(iz)
    IF(iz .gt. 0) dslk  = AMIN1(dsl_z(iz),vk*ABS(z(iz))/csmag)
    ! NO STABILITY CORRECTED LENGTH SCALES
    DO j=iys,iye
      DO i=1,nnx
        alk(i,j,iz) = dslk
      END DO
    END DO
    IF (dzw_i(izp1) == 0.0 .OR. dzu_i(izp1) == 0.0) THEN
        print*, 'Zero denominator at iz=', iz, ' -> dzw_i(izp1)=', dzw_i(izp1), ' dzu_i(izp1)=', dzu_i(izp1)
    ENDIF
    ! GET STRAINS and compute viscosity
    DO j=iys,iye
      DO i=1,nnx
        s11 = (ux(i,j,iz))**2
        s22 = (vy(i,j,iz))**2
        wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
        s33 =  (wzp)**2
        s12 =  (0.5 * (uy(i,j,iz) + vx(i,j,iz)))**2
        s13 =  (0.5 * ((((u(i,j,izp1) - u(i,j,iz)) * dzu_i(izp1) + wx(i,j,iz)))))**2
        s23 =  (0.5 * ((v(i,j,izp1) - v(i,j,iz)) * dzu_i(izp1) + wy(i,j,iz)))**2
        print *, "s11", s11, "s22", s22, "s33", s33, "s12", s12, "s13", s13, "s23", s23
        sij2(i,j,iz) = s11 + s22 + s33 + 2 * s12 + 2 * s13 + 2 * s23
        vis_m(i,j,iz)  = (csmag*dslk)**2 * SQRT(2 * sij2(i, j, iz))
        vis_s(i,j,iz)  = vis_m(i,j,iz)
        vis_sv(i,j,iz) = vis_s(i,j,iz)
      ENDDO
    ENDDO
    print*, "iz", iz, "vis_m", vis_m(1,1,iz), "vis_s", vis_s(1,1,iz), "vis_sv", vis_sv(1,1,iz), "dslk", dslk, "sij2", sij2(1,1,iz)
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