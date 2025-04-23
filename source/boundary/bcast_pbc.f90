SUBROUTINE bcast_pbc
! SEND UPPER BC TO OTHER PROCESSORS FOR FFT SOLUTION OF PRESSURE

  USE pars
  USE fields

  INCLUDE 'mpif.h'

  INTEGER :: istatus(mpi_status_size)
  INTEGER :: irow_r

  IF(numprocs /= 1) THEN
    irow_r = MOD(myid,ncpu_s)
    irow_t = is_s(numprocs-1) + irow_r
    num = nnx*(iye+1-iys)

    ! CHECK WHICH ROW MYID IS IN
    IF(iss /= is_s(numprocs-1)) THEN
      ! NOT IN TOP ROW, RECEIVE FROM TOP
      CALL mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,mpi_comm_world,     &
      istatus,ierr)
    ELSE
      ! MYID IS IN THE TOP ROW, SEND BELOW
      DO l=irow_r,irow_t-ncpu_s,ncpu_s
        CALL mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,mpi_comm_world,ierr)
      END DO
    ENDIF

    ! REPEAT FOR NEXT VARIABLE
    IF(iss /= is_s(numprocs-1)) THEN
      ! NOT IN TOP ROW, RECEIVE FROM TOP
      CALL mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,mpi_comm_world,    &
      istatus,ierr)
    ELSE
      ! MYID IS IN THE TOP ROW, SEND BELOW
      DO l=irow_r,irow_t-ncpu_s,ncpu_s
        CALL mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,mpi_comm_world,ierr)
      END DO
    ENDIF
  ELSE
    CONTINUE
  ENDIF

  RETURN
END SUBROUTINE
