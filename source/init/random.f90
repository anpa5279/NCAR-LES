SUBROUTINE random
! GEOSTROPHIC WINDS DESIGNED FOR COMPARISON CASE

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  REAL :: psi(nnx,iys:iye), psix(nnx,iys:iye), psiy(nnx,iys:iye,izs:izs),   &
        uxx(nnx,iys:iye), vyy(nnx,iys:iye,izs:izs)

  ! NOTE SET NMATCH IN SR. ISO SO THAT IT IS COMPATIBLE WITH CONDITIONS HERE
  DO iz=1,nnz
    ug(iz)   = ugcont
    vg(iz)   = vgcont
    divz(iz) = 0.0
  ENDDO

  izi = (50*nnz)/100
  zi  = z(izi)

  z_lower = zi - 50.0
  t_lower = 300.0
  z_upper = zi + 50.0
  t_upper = 308.0
  slope   = (t_upper - t_lower)/(z_upper - z_lower)

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz) = ugcont-ugal
        v(ix,iy,iz) = vgcont
        w(ix,iy,iz) = 0.0
        e(ix,iy,iz) = 0.0
      ENDDO
    ENDDO

    IF(z(iz) <= z_lower) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,1,iz) = t_lower
        ENDDO
      ENDDO
    ELSEIF(z(iz) >= z_upper) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,1,iz) = t_upper + (zz(iz) - z_upper)*dtdzf(1)
        ENDDO
      ENDDO
    ELSE
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,1,iz) = t_lower + slope*(zz(iz) - z_lower)
        ENDDO
      ENDDO
    ENDIF

    DO iy=iys,iye
      DO ix=1,nnx
        w(ix,iy,iz)   = 0.
        r1(ix,iy,iz)  = 0.
        r2(ix,iy,iz)  = 0.
        r3(ix,iy,iz)  = 0.
        r4(ix,iy,1,iz)= 0.
        r5(ix,iy,iz)  = 0.
      ENDDO
    ENDDO
  ENDDO

  ! SET INITIAL RANDOM FIELD TO BE DIVERGENCE FREE
  idum = -1 - myid
  DO iz=izs,ize

    ! AMPV AND AMPT ARE MAX AMPLITUDES OF RANDOM VELOCITY AND TEMPERATURE FIELDS
    ! MAKE SURE AMPV IS SET IF FREE CONVECTION SO THAT WE HAVE MOTIONS AT FIRST
    ! TIME STEP
    ampv = 0.0
    ampv = 0.001
    ampt = 0.10

    ! SIMPLE RANDOM FIELD SCALED BETWEEN -0.5 AND 0.5
    sum_psi = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        psi(ix,iy) = ran1(idum)
        sum_psi = sum_psi + psi(ix,iy)
      ENDDO
    ENDDO

    sum_psi = sum_psi*fnxy
    CALL mpi_sum_xy(sum_psi,myid,iss,ise,1)

    DO iy=iys,iye
      DO ix=1,nnx
        psi(ix,iy) = psi(ix,iy) - sum_psi
        psix(ix,iy)     = psi(ix,iy)
        psiy(ix,iy,izs) = psi(ix,iy)
      ENDDO
    ENDDO

    CALL xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
    CALL yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e, &
          iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)

    vmaxx = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        vmag = SQRT(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
        IF(vmag > vmaxx) vmaxx = vmag
      ENDDO
    ENDDO

    facv = ampv/vmaxx
    IF (z(iz) <= 50.0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
          v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
          t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
        ENDDO
      ENDDO
    ENDIF

    IF(z(iz) <= 250.0) THEN
      DO iy=iys,iye
        DO ix=1,nnx
          e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
        ENDDO
      ENDDO
    ENDIF

    ! CHECK DIVERGENCE OF INITIAL FIELD
    DO iy=iys,iye
      DO ix=1,nnx
        uxx(ix,iy) = u(ix,iy,iz)
        vyy(ix,iy,izs) = v(ix,iy,iz)
      ENDDO
    ENDDO

    CALL xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
    CALL yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),nnx,nny,ixs,ixe,ix_s,ix_e,  &
          iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)

    DO iy=iys,iye
      DO ix=1,nnx
        divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
      ENDDO
    ENDDO
    divz(iz) = divz(iz)*fnxy
  ENDDO

  CALL mpi_sum_z(divz(1),i_root,myid,nnz,1)

  WRITE(nprt,6000)
  WRITE(nprt,6100) (iz,divz(iz),iz=izs,ize)

  ! FIX FOR BAROCLINIC AND SUBSIDENCE EFFECTS
  RETURN

! FORMAT
6000  FORMAT(' check of divergence for initial state',/,' iz ',5x,' divergence')
6100  FORMAT(i5,e15.6)

END SUBROUTINE
