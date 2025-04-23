SUBROUTINE restart
! GET RESTART FILE FROM LOCAL DIRECTORY

  USE pars
  USE fields
  USE con_data
  USE con_stats

  CHARACTER*80 path_res_c
  LOGICAL :: there

  ! CHECK IF FILE IS THERE
  INQUIRE(file=path_res,exist=there)

  IF(there) THEN
    IF(l_root) WRITE(6,6001) path_res
  ELSE
    IF(l_root) WRITE(6,6005) path_res
    STOP
  ENDIF

  ! GET CONSTANT FILE
  iloc = INDEX(path_res,' ')
  path_res_c = path_res(1:iloc-1)//'.con'

  INQUIRE(file=path_res_c,exist=there)

  IF(there) THEN
    IF(l_root) WRITE(6,6002) path_res_c
  ELSE
    IF(l_root) WRITE(6,6006) path_res_c
    STOP
  ENDIF

  OPEN(nvelc,err=200,file=path_res_c,form='unformatted',status='old')

  CALL read_res

  RETURN

  ! PROCESS ERRORS
  200 CONTINUE
  WRITE(6,9001) path_res_c, nvelc
  CALL mpi_finalize(ierr)
  STOP

! FORMAT
6001  FORMAT(' SR. RESTART: FILE READ = ',A80)
6002  FORMAT(' SR. RESTART: CONSTANT FILE READ = ',A80)
6005  FORMAT(' 6005, SR. RESTART: cannot find restart file = ',a80)
6006  FORMAT(' 6005, SR. RESTART: cannot find constant file = ',a80)
9000  FORMAT(' 9000, SR. RESTART: cannot open file =',a80,/,                &
            ' to unit number = ',i2)
9001  FORMAT(' 9001, SR. RESTART: cannot open file =',a80,/,                &
            ' to unit number = ',i2)

END SUBROUTINE
