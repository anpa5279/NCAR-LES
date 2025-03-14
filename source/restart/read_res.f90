SUBROUTINE read_res
! READ RESTART FILE INCLUDING CONSTANT FILE
! CHANGED FOR IYS:IYE

  USE pars
  USE fields
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'

  INTEGER :: status(mpi_status_size), ierr
  INTEGER(kind=mpi_offset_kind) :: offset, disp
  INTEGER(kind=k8)              :: nsize, nsize2
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp

  ALLOCATE(temp(nvar,nnx,iys:iye))

  ! OPEN FILE
  CALL mpi_file_open(mpi_comm_world, path_res,                              &
        mpi_mode_create+mpi_mode_rdwr,mpi_info_null, nvel, ierr)

  ! SET FILE VIEW
  disp = 0
  CALL mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,'native',            &
        mpi_info_null,ierr)

  ! READ 3D FIELDS
  nsize  = INT(nvar,k8)*nnx*nny
  nsize2 = INT(nvar,k8)*nnx*(iys-1)
  n_read = nvar*nnx*(iye+1-iys)

  DO k=izs,ize
    offset = INT((k-1),k8)*nsize + nsize2
    CALL mpi_file_read_at_all(nvel,offset,temp,n_read,mpi_real8,status,ierr)

    IF (ierr /= 0) goto 9992

    DO j=iys,iye
      DO i=1,nnx
        u(i,j,k) = temp(1,i,j)
        v(i,j,k) = temp(2,i,j)
        w(i,j,k) = temp(3,i,j)
        e(i,j,k) = temp(nvar,i,j)
      ENDDO
    ENDDO

    DO is = 1,nscl
      DO j = iys,iye
        DO i = 1,nnx
          t(i,j,is,k) = temp(3+is,i,j)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! CLOSE FILE
  CALL mpi_file_close(nvel, ierr)

  DEALLOCATE(temp)

  ! EVERY MPI PROCESS READS CONSTANT FILE
  REWIND(nvelc)
  READ(nvelc,err=9993) c_c, c_hurr, c_s, case

  IF(iti <= 1)THEN
    time = 0.1
    dt = 0.1 !initial time step size
  ENDIF

  CLOSE(nvelc)

  IF(l_root) WRITE(6,4001) case

  ! SPECIAL RESTART CONDITIONS
  ! SET CASE NAME TO CASE INPUT
  case   = case_inp
  IF(l_root) WRITE(6,4002) case_inp, utau, utausv

  ! IF NEW VIS MODEL SET MATCH POINT FOR OUTER GRID
  nmatch = nnz
  utau = utausv

  ! REDEFINE CASE ID TO INPUT VALUE
  IF(l_root) WRITE(6,4012) time
  IF(l_root) WRITE(6,4013) qstar(1), nmatch, case

  CALL get_dz

  RETURN

  9992 CONTINUE
  WRITE(6,6100) nvel,iz

  CALL mpi_finalize(ierr)
  STOP

  9993 CONTINUE
  WRITE(6,6200) nvelc

  CALL mpi_finalize(ierr)
  STOP

! FORMAT
4001  FORMAT(' 4001, SR. RESTART: case from restart = ',a3)
4002  FORMAT(' 4002, SR. RESTART:',/,                                       &
            ' files will be saved with case name = ',a3,/,' utau = ',e15.6, &
            ' utausv = ',e15.6)
6100  FORMAT(' SR. READ_RES: error reading file on unit number = ',i2,/,    &
            '               at iz = ',i4)
6200  FORMAT(' SR. READ_RES:',/,                                            &
            '    error reading constant file on unit number = ',i2)
4012  FORMAT(' SR. RESTART: restart completed at T=',e15.6)
4013  FORMAT('    after restart qstar = ',e15.6,' nmatch = ',i5,' case = ',a3)

END SUBROUTINE
