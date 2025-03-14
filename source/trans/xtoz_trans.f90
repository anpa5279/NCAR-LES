SUBROUTINE xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,iys,iye,izs,ize,iz_s,     &
  iz_e,myid,ncpu_s,numprocs)
! only called in solve/pressure
! TRANSPOSE ARRAY F(NX,IYS:IYE,IZS-1:IZE-1) --> G(0:NZ+1,IYS:IYE,IXS:IXE)

  INCLUDE 'mpif.h'

  INTEGER :: istatus(mpi_status_size)
  REAL :: f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
  REAL :: ft(nx*(iye+1-iys)*(ize-izs+1)), gt(nz*(ixe+1-ixs)*(iye-iys+1))
  INTEGER :: ix_s(0:numprocs-1), ix_e(0:numprocs-1), iz_s(0:numprocs-1),    &
        iz_e(0:numprocs-1)

  jk = (ize - izs + 1)*(iye - iys + 1)
  ij = (ixe - ixs + 1)*(iye - iys + 1)

  ! GET STARTING LOCATION
  iss = myid - (myid/ncpu_s)*ncpu_s

  DO i=iss,numprocs-1,ncpu_s
    nsend = (ix_e(i) - ix_s(i) + 1)*jk
    nrecv = (iz_e(i) - iz_s(i) + 1)*ij
    IF(i == myid) THEN
      CALL send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),iys,iye,iz_s(myid),iz_e(myid))
    ELSE
      CALL send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),iys,iye,iz_s(myid),iz_e(myid))
      CALL mpi_sendrecv(ft(1),nsend,mpi_real8,i,1,gt(1),nrecv,mpi_real8,i,  &
            1,mpi_comm_world,istatus,ierr)
    ENDIF
    CALL recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),iys,iye,iz_s(i),iz_e(i))
  ENDDO

  RETURN
END
