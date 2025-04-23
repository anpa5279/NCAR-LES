SUBROUTINE pbltop(itop)
! ---------- get estimate of pbl top
!            method = 0, min of wt flux
!                        (good for buoyancy cases)
!            method = 1, uw flux less than critical value
!                        (good for ekman cases)
!            method = 2, running t average exceeds criterion
!                        (good for neutral cases with capping
!                         inversions)
!            method = 3, MAXimum gradient in temperature field
!                        (good for finding local zi see jas paper)
!                        with minimum search height (sr. setup)
! ------------ IF method USEs average statistics THEN only root
!              process need find zi

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: trun(maxnz)
  INCLUDE 'mpif.h'
  REAL :: gradloc(2,nnx,nny), gradmax(2,nnx,nny)
  EXTERNAL get_zi

  IF(method <= 2 .AND. l_root) THEN
    sgn = 1.0
    IF (method <= 0 .OR. method > 2) THEN
      itop=1
      wttot=wtle(1,1)+wtsb(1,1)
      wtmin=wttot*sgn
      DO iz=2,nnz
        wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
        IF (wttot<=wtmin) THEN
          itop=iz
          wtmin=wttot
        ENDIF
      ENDDO
      zi=z(itop)
    ELSE IF (method == 1) THEN
      itop = 1
      crit = 0.05
      uwsf = utau*utau
      DO iz=1,nnzm1
        uwtot = (uwle(iz) + uwsb(iz))**2 + (vwle(iz) + vwsb(iz))**2
        uwtot = SQRT(uwtot)
        IF(uwtot/uwsf > crit) THEN
          itop=iz
        ENDIF
      ENDDO
      zi=z(itop)
    ELSE IF (method == 2) THEN
      trun(1) = txym(1,1)
      DO iz=2,nnz
        weight = z(iz-1)/z(iz)
        trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
      ENDDO
      itop = 1
      tcrit = 0.1
      DO iz=2,nnz
        IF(txym(iz,1) > (trun(iz) + tcrit)) THEN
          itop = iz
          EXIT
        ENDIF
        EXIT
      ENDDO
      zi=z(itop)
    ENDIF

    DO iy=1,nny
      DO ix=1,nnx
        gradmax(2,ix,iy) = zi
      ENDDO
    ENDDO

  ! USE GRADIENT METHOD, EVERY PROCESS COMPUTES
  ELSEIF(method == 3) THEN

    ! SIMILAR TO ZEROING IN STAT ARRAY IN SR. MEAN_STAT
    DO iy=1,nny
      DO ix=1,nnx
        gradloc(1,ix,iy) = 0.0
        gradloc(2,ix,iy) = z(iz_min)
      ENDDO
    ENDDO

    ! NOW ALL Z IN THIS PROCESS
    IF(iz_min <= ize) THEN
      DO iz=MAX(izs,iz_min),ize
        izp1 = iz + 1
        DO iy=iys,iye
          DO ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            IF(grad > gradloc(1,ix,iy)) THEN
              gradloc(1,ix,iy) = grad
              gradloc(2,ix,iy) = z(iz)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! ALTERNATE VERSION USING ALREADY DEFINED FUNCTION IN MPI
    ! PASSES 2 REAL8 VARIABLES
    CALL mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2DOuble_precision,          &
          mpi_maxloc,i_root,mpi_comm_world,ierror)

    ! GET AVERAGE ON ROOT PROCESSES
    IF(l_root) THEN
      zi_avg = 0.0
      DO iy=1,nny
        DO ix=1,nnx
          zi_avg = zi_avg + gradmax(2,ix,iy)
        ENDDO
      ENDDO
      zi = zi_avg*fnxy
    ENDIF
  ENDIF

  ! SEND AVERAGE ZI EVERYWHERE
  CALL mpi_bcast(zi,1,mpi_real8,i_root,mpi_comm_world,ierr)

  DO iz=1,nnz
    IF(zi <= z(iz) .AND. zi > z(iz+1)) itop = iz
  ENDDO

  RETURN

! FORMAT
7001  FORMAT(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,' itop = ',i3)

END SUBROUTINE
