SUBROUTINE xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,iz1,    &
  iz2,myid,ncpu_s,np)
!called in deriv/yd_mpi and fft/fft2d_mpi
!TRANSPOSE ARRAY F(NX,IYS:IYE,IZ1:IZ2) --> G(NY,IXS:IXE,IZ1:IZ2)

  INCLUDE 'mpif.h'

  INTEGER :: istatus(mpi_status_size)
  REAL :: f(nx,iys:iye,iz1:iz2),g(ny,ixs:ixe,iz1:iz2)
  REAL :: ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
  INTEGER :: ix_s(0:np-1), ix_e(0:np-1),iy_s(0:np-1), iy_e(0:np-1)

  jk = (iye - iys + 1)*(iz2 - iz1 + 1)
  ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

  !GET CPUS ON SLAB FOR MYID
  islab = myid/ncpu_s
  iss   = islab*ncpu_s
  ise   = iss + ncpu_s - 1

  DO i=iss,ise
    nsend = (ix_e(i) - ix_s(i) + 1)*jk
    nrecv = (iy_e(i) - iy_s(i) + 1)*ik
    IF(i == myid) THEN
      CALL send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),iy_s(myid),iy_e(myid),iz1,iz2)
    ELSE
      CALL send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),iy_s(myid),iy_e(myid),iz1,iz2)
      CALL mpi_sendrecv(ft(1),nsend,mpi_real8,i,1,gt(1),nrecv,mpi_real8,i,  &
            1,mpi_comm_world,istatus,ierr)
    ENDIF

    CALL recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),iy_s(i),iy_e(i),iz1,iz2)
  ENDDO

  RETURN
END SUBROUTINE
