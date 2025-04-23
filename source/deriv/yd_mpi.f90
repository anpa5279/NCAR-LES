SUBROUTINE yd_mpi(ay,trigy,yk,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,    &
  iz1,iz2,myid,ncpu,np)
! GET MULTIPLE Y DERIVATIVES USING FFTPACK ROUTINES AND MPI
! USE FFTPACK STORAGE A_O, (A1,B1), (A2,B2), ...
! ASSUMES THAT FIRST WAVENUMBER YK(1) = 0.0
! WAVE NUMBERS ARE NORMALIZED BY NUMBER OF POINTS, NY

  REAL :: yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
  REAL :: ayt(ny,ixs:ixe,iz1:iz2)
  INTEGER :: ix_s(0:np-1), ix_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)

  CALL xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,iz1,iz2, &
        myid,ncpu,np)

  DO iz=iz1,iz2
    DO ix=ixs,ixe
      CALL rfftf(ny,ayt(1,ix,iz),trigy)

      ii = 1
      ayt(1,ix,iz)  = 0.0
      ayt(ny,ix,iz) = 0.0
      DO iy=2,ny-1,2
        ii              = ii + 1
        temp            = ayt(iy,ix,iz)
        ayt(iy,ix,iz)   = -yk(ii)*ayt(iy+1,ix,iz)
        ayt(iy+1,ix,iz) = yk(ii)*temp
      ENDDO

      CALL rfftb(ny,ayt(1,ix,iz),trigy)
    ENDDO
  ENDDO

  CALL ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,iz1,iz2, &
        myid,ncpu,np)

  RETURN
END SUBROUTINE yd_mpi
