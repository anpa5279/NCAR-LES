SUBROUTINE dear_vis(alk)
! GET DEARDORFF STABILITY CORRECTED LENGTH SCALES, POSSIBLE NEW VIS MODEL

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: alk(nnx,iys:iye,izs-1:ize+1)

  DO iz=izs-1,ize+1
    izp1 = iz + 1
    dslk  = dsl_z(iz)
    IF(iz > 0) dslk  = AMIN1(dsl_z(iz),vk*ABS(z(iz))/csmag)
    almin = almin_c*dsl_z(iz)
    IF(iz == 0 .or. iz == nnz+1) THEN
      dfack = 1.0
    ELSE
      dfack = dfac(iz)
    ENDIF

    IF(ivis == 1 .AND. iz <= nmatch) THEN
      ! NO STABILITY CORRECTED LENGTH SCALES
      DO j=iys,iye
        DO i=1,nnx
          alk(i,j,iz) = dslk
        END DO
      END DO
    ELSE
      DO j=iys,iye
        DO i=1,nnx
          alk(i,j,iz) = dslk
          stab = batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1)
          IF(stab>stabmin) THEN
            als = stab_c*SQRT(e(i,j,iz)/stab)
            alk(i,j,iz) = AMIN1(dslk,als)
          ENDIF
          alk(i,j,iz)  = AMAX1(almin,alk(i,j,iz))
        ENDDO
      ENDDO
    ENDIF

    DO j=iys,iye
      DO i=1,nnx
        vis_m(i,j,iz)  = ck*alk(i,j,iz)*SQRT(e(i,j,iz))*dfack
        vis_s(i,j,iz)  = (1.+2.*alk(i,j,iz)/dslk)*vis_m(i,j,iz)
        vis_sv(i,j,iz) = vis_s(i,j,iz)
      ENDDO
    ENDDO

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
