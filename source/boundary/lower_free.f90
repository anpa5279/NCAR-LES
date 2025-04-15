SUBROUTINE lower_free(it)
!!!CHECK LOWER FREE FLAG (if parameter ifree=1) in les_mpi.f90. if ifree==0 and iss==0, then lower.f90 instead of this
!SETUP LOWER BC FOR FREE CONVECTION WHERE EACH PROCESSOR APPLIES LOG-LAW AT
!SEVERAL (IX,IY) FOR IZ = 1
!INDEX F(.,.,2) INDICATES LOWER

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'

  REAL :: u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
  REAL :: sbuf(2+2*nscl,mxs:mxe,iys:iye)
  REAL :: rbuf((2+2*nscl)*nnx*(iye+1-iys))

  !BROADCAST LEVEL 1 DATA EVERYWHERE
  IF(iss == 0) THEN
    DO iy=iys,iye
      DO ix=1,nnx
        u_level1(ix,iy,1) = u(ix,iy,1)
        u_level1(ix,iy,2) = v(ix,iy,1)
      ENDDO
    ENDDO

    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  num = nnx*(iye + 1 - iys)*(2+nscl)

  !SEND ALL OF ROOT DATA TO OTHER PROCESSORS
  CALL mpi_send_root(u_level1(1,iys,1),num,myid,numprocs,ncpu_s)

  !EVERY TASK GETS THEIR OWN FLUXES AND SURFACE SCALARS
  CALL suft2(u_level1,it)

  !SEND SURFACE SCALARS AND MOMENTUM FLUXES BACK TO ROOT(S)
  IF(numprocs /= 1) THEN !numprocs= total number of processes OR size of cluster
    DO iy=iys,iye
      DO ix=mxs,mxe
        sbuf(1,ix,iy)  = tau13m(ix,iy)
        sbuf(2,ix,iy)  = tau23m(ix,iy)
      ENDDO
    ENDDO

    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=mxs,mxe
          sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
          sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
        ENDDO
      ENDDO
    ENDDO

    irow_r = MOD(myid,ncpu_s)
    IF(myid >= ncpu_s) THEN !if the process rank is greater than or equal to ncpu_s (8)
      num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
      CALL mpi_send(sbuf(1,mxs,iys),num,mpi_REAL8,irow_r,1,mpi_comm_world,ierr)
    ELSE
      DO l=irow_r+ncpu_s,numprocs-1,ncpu_s
        num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
        CALL mpi_recv(rbuf(1),num,mpi_REAL8,l,1,mpi_comm_world,istatus,ierr)
        CALL f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,tau13m,tau23m,   &
              taut3m,t_grnd)
      ENDDO
    ENDIF

  ELSEIF(numprocs == 1) THEN !numprocs= total number of processes OR size of cluster
    CONTINUE
  ENDIF

!ONLY FOR ROOT ROW = 0
!GET SUMS OF SURFACE CONDITIONS AND SET SURFACE BOUNDARY CONDITIONS
  IF(iss == 0) THEN
    buf(1) = 0.0
    buf(2) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        buf(1) = buf(1) + tau13m(ix,iy)
        buf(2) = buf(2) + tau23m(ix,iy)
      ENDDO
    ENDDO

    DO iscl=1,nscl
      buf(2+iscl)      = 0.
      buf(2+nscl+iscl) = 0.
      DO iy=iys,iye
        DO ix=1,nnx
          buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
          buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
        ENDDO
      ENDDO
    ENDDO

    CALL mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)

    uwsfc = buf(1)*fnxy
    vwsfc = buf(2)*fnxy
    DO iscl=1,nscl
      wtsfc(iscl) = buf(2+iscl)*fnxy
      tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
    ENDDO

    iz   = 1
    izm1 = iz - 1
    dz_i = dzu_i(iz) !inverse dz

    DO iy=iys,iye
      DO ix=1,nnx
        ebc(ix,iy,2)=AMAX1(e(ix,iy,iz),sml_eg)
        wbc(ix,iy,2)= 0.0
        pbc(ix,iy,2) = 0.0
        pbc2(ix,iy,2) = 0.0
      ENDDO
    ENDDO

    DO iy=iys,iye
      DO ix=1,nnx
        dudz     = 2.*u(ix,iy,iz)*dz_i
        dvdz     = 2.*v(ix,iy,iz)*dz_i
        ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
        vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
      ENDDO
    ENDDO

    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
          tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
        ENDDO
      ENDDO
    ENDDO

    !INITIALIZE U,V,W,T AND DERIVATIVES AT IZM1
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,izm1)  = ubc(ix,iy,2)
        v(ix,iy,izm1)  = vbc(ix,iy,2)
        w(ix,iy,izm1)  = wbc(ix,iy,2)
        r3(ix,iy,izm1) =  0.0
        e(ix,iy,izm1)  = ebc(ix,iy,2)
        ux(ix,iy,izm1) = 0.0
        uy(ix,iy,izm1) = 0.0
        vx(ix,iy,izm1) = 0.0
        vy(ix,iy,izm1) = 0.0
        wx(ix,iy,izm1) = wbc(ix,iy,2)
        wy(ix,iy,izm1) = wbc(ix,iy,2)
      ENDDO
    ENDDO

    DO iscl=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  RETURN
END SUBROUTINE
