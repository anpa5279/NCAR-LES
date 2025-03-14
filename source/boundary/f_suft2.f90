SUBROUTINE f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,tau13m,tau23m,taut3m,t_grnd)
  !only used for boundary/lower_free.f90. only used if ifree=1 and if numprocs does not equal 1
  IMPLICIT NONE

  !FILL SURFACE ARRAYS ON ROOT PROCESSORS
  INTEGER :: iscl , ix , iy , iye , iys , mxe , mxs , nnx , nscl
  !nscl= number of scalars and vars (nscl=1 is only physics)
  REAL :: rbuf(2+2*nscl,mxs:mxe,iys:iye) !how rbuf is a 3D matrix? rows= scalars/ tracking chemical quantities, cols= x, cols2= y (remember parallelizing in the z)
  REAL ::                                                                   &
          tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),       & !2D matrices
          taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl) !3D matrices
          !tau = shear stress
  DO iy=iys,iye !y
    DO ix=mxs,mxe !x
      tau13m(ix,iy) = rbuf(1,ix,iy)
      tau23m(ix,iy) = rbuf(2,ix,iy)
    ENDDO
  ENDDO

  DO iscl=1,nscl !scalars
    DO iy=iys,iye !y
      DO ix=mxs,mxe !x
        taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
        t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
