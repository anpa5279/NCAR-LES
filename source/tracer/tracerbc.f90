module tracerbc
  use con_data, only: dz,dx,dy
  use con_stats, only: z, zz
  use pars, only: iys,iye,izs,ize,izi,nnx,nnz,nny,nscl
  use fields, only: t
  implicit none
  integer, parameter :: flg_debug = 0
  integer, dimension(nscl) :: ictype
  integer, dimension(2,nscl) :: bnd
  integer, dimension(2) :: bnds
  integer, dimension(3,nscl) :: point
  integer, dimension(3) :: points
  real, dimension(nscl) :: val
  contains

! iscl      : scalar number (temperature is always iscl=1)
! tau       : reaction time scale
! ictype    : initial condition (0 = air-sea flux, 1 = horiz. band, 
!                                2 = vertical band in x, 3 = vertical band in y,
!                                4 = point source, 5 = vertical gradient, 
!                                6 = horiz. gradient in x, 7 = horiz. gradient in y)
! kconst    : air-sea interface flux rate (if 0 then there is no flux across top boundary)
! val       : value of initial finite or source band/point
! np        : width of initial finite or source band
! zt        : upper/left most level or finite or source band
! bnd       : 
! point     : 
! rmodel    : reaction model type (0 = no reaction, 1 = single tracer decay/growth,
!                                  2 = two tracers decay/growth, 3 = carbonate chemistry)
! rdorg     : reaction decay or growth (0 = decaying tracer, 1 = growing tracer)
! rpartner  : reaction partner (iscl number of coupled tracer for reaction, 
!                               0 = no coupled tracer)

    subroutine applytracerbc(it)
      integer, intent(in) :: it
      integer :: iscl, np, zt
      real :: ta, vals

      !! tracer 1 - temperature
      iscl = 1;  
      
      !! tracers 
      iscl = 2;   
      ictype(iscl) = 1;   val(iscl) = 14.0
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 3;   
      ictype(iscl) = 1;   val(iscl) = 0.18
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 4;   
      ictype(iscl) = 1;   val(iscl) = 0.01
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 5;   
      ictype(iscl) = 1;   val(iscl) = 0.12
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 6;   
      ictype(iscl) = 1;   val(iscl) = 2.0
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 7;   
      ictype(iscl) = 1;   val(iscl) = 3.0e-3
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 8;  
      ictype(iscl) = 1;   val(iscl) = 0.025
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 9;   
      ictype(iscl) = 1;   val(iscl) = 180.0
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 10;   
      ictype(iscl) = 1;   val(iscl) = 0.031
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 11;   
      ictype(iscl) = 1;   val(iscl) = 1.5e-3
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 12;   
      ictype(iscl) = 1;   val(iscl) = 0.28
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 13;   
      ictype(iscl) = 1;   val(iscl) = 3.5e-3
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 14;   
      ictype(iscl) = 1;   val(iscl) = 1.5e-4
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 15;  
      ictype(iscl) = 1;   val(iscl) = 230.0
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);
      
      iscl = 16;   
      ictype(iscl) = 1;   val(iscl) = 0.08
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 17;   
      ictype(iscl) = 1;   val(iscl) = 1.5
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      iscl = 18;  
      ictype(iscl) = 1;   val(iscl) = 0.18
      np = nnz+2;  zt = 0; bnd(:,iscl) = znptobnd(zt,np);

      do iscl = 2,nscl
         bnds=bnd(:,iscl); vals=val(iscl); points=point(:,iscl);
         if (it.eq.1) then
            if (ictype(iscl).eq.1) call hbndsource(iscl,bnds,vals);
            !if (ictype(iscl).eq.2) call vxbndsource(iscl,bnds,vals);
            !if (ictype(iscl).eq.3) call vybndsource(iscl,bnds,vals);
            !if (ictype(iscl).eq.5) call vgradsource(iscl,bnds,vals);
            !if (ictype(iscl).eq.6) call hxgradsource(iscl,bnds,vals);
            !if (ictype(iscl).eq.7) call hygradsource(iscl,bnds,vals);
!            call bfm_profiles
         endif
         !if (ictype(iscl).eq.4) call pointsource(iscl,points,vals); 
         bnds = 0; vals = 0; points = 0;
      enddo
      
      if(flg_debug == 1) then
          open(13, file='tracerbc.txt',access='append')
          write(13,'(A)') '------------------------'
          write(13,'(A,i3)') 'RUNNING FOR IT= ',it
          write(13,'(A,f9.6)') 'Z for 5m above H', z(izi)+5.0
          write(13,'(A,f9.6)') 'Z for 5m below H', z(izi)-5.0
          close(13)
      end if

    end subroutine

    subroutine bfm_profiles
      integer :: ix,iy,iz,iscl
      do iz=izs,ize
         do iy=iys,iye
            do ix=1,nnx
!               do iscl=2,nscl
!                  t(ix,iy,iscl,iz)=1.0
!               end do
               t(ix,iy,2,iz) = -7.0*tanh(1.25*abs(z(iz))-100.0)+7.0
               t(ix,iy,3,iz) = -0.09*tanh(1.25*abs(z(iz))-100.0)+0.09
               t(ix,iy,4,iz) = -0.005*tanh(1.25*abs(z(iz))-100.0)+0.005
               t(ix,iy,5,iz) = -0.06*tanh(1.25*abs(z(iz))-100.0)+0.06

               t(ix,iy,6,iz) = -1.0*tanh(1.25*abs(z(iz))-100.0)+1.0 +1.0e-10
               t(ix,iy,7,iz) = -(1.5e-3)*tanh(1.25*abs(z(iz))-100.0)+1.5e-3
               t(ix,iy,8,iz) = -0.0125*tanh(1.25*abs(z(iz))-100.0)+0.0125

               t(ix,iy,9,iz) = 180.0*exp(-abs(zz(iz)))**0.03 + 90.0
!               if(abs(z(iz)).le.100.0)then
                  t(ix,iy,10,iz) = -0.015*tanh(1.25*abs(z(iz))-100.0)+0.015
                  t(ix,iy,11,iz) = -(0.75e-3)*tanh(1.25*abs(z(iz))-100.0) + 0.75e-3
                  t(ix,iy,12,iz) = -(0.14)*tanh(1.25*abs(z(iz))-100.0) + 0.14
                  t(ix,iy,13,iz) = -(1.75e-3)*tanh(1.25*abs(z(iz))-100.0)+1.75e-3
                  t(ix,iy,14,iz) = -(0.75e-4)*tanh(1.25*abs(z(iz))-100.0)+0.75e-4
!                  t(ix,iy,16,iz) = 0.0
!                  t(ix,iy,17,iz) = 0.0
!               else
!                  t(ix,iy,10,iz) = 0.0
!                  t(ix,iy,11,iz) = 0.0
!                  t(ix,iy,12,iz) = 0.0
!                  t(ix,iy,13,iz) = 0.0
!                  t(ix,iy,14,iz) = 0.0
                  t(ix,iy,16,iz) =0.02*tanh(1.25*abs(z(iz))-100.0)+0.0242
                  t(ix,iy,17,iz) =0.5*tanh(1.25*abs(z(iz))-100.0)+0.86
!               endif

               t(ix,iy,15,iz) = 20.0*exp(-abs(z(iz)))**0.03+194.0
               t(ix,iy,18,iz) = -0.09*tanh(1.25*z(iz)-100.0)+0.09
               
            end do
         end do
      end do
    end subroutine

    subroutine hbndsource(iscl, bnds, vals)
      integer, intent(in) :: iscl
      integer, intent(in), dimension(2) :: bnds
      real, intent(in) :: vals
      integer :: ix,iy,iz
      do iz=bnds(1),bnds(2)
         do iy=iys,iye
            do ix=1,nnx
               if ((iz >= izs) .and. (iz <= ize)) then
                     t(ix,iy,iscl,iz) = vals
               endif
            end do
         end do
      end do

    end subroutine

!!$    subroutine vxbndsource(iscl, bnds, vals)
!!$      integer, intent(in) :: iscl
!!$      integer, intent(in), dimension(2) :: bnds
!!$      integer, intent(in) :: vals
!!$      integer :: ix,iy,iz
!!$      do iy=iys,iye
!!$         do iz=izs,ize
!!$            do ix=bnds(1),bnds(2)
!!$               t(ix,iy,iscl,iz) = vals
!!$            end do
!!$         end do
!!$      end do
!!$
!!$    end subroutine

!!$    subroutine vybndsource(iscl, bnds, vals)
!!$      integer, intent(in) :: iscl
!!$      integer, intent(in), dimension(2) :: bnds
!!$      integer, intent(in) :: vals
!!$      integer :: ix,iy,iz
!!$      do iy=bnds(1),bnds(2)
!!$         do iz=izs,ize
!!$            do ix=1,nnx
!!$               if ((iy >= iys) .and. (iy <= iye)) then
!!$                     t(ix,iy,iscl,iz) = vals
!!$               endif
!!$            end do
!!$         end do
!!$      end do
!!$
!!$    end subroutine

!    subroutine vgradsource(iscl, bnds, vals)
!      integer, intent(in) :: iscl
!      integer, intent(in), dimension(2) :: bnds
!      real, intent(in) :: vals
!      integer :: ix,iy,iz,zi

!      zi  = z(bnds(2))
!      do iy=iys,iye
!         do iz=bnds(1),bnds(2)
!            do ix=1,nnx
!               if ((iz >= izs) .and. (iz <= ize)) then
!                  t(ix,iy,iscl,iz) = (vals/zi)*(zi-zz(iz))
!               endif
!            end do
!         end do
!      end do

!    end subroutine

!!$    subroutine hxgradsource(iscl, bnds, vals)
!!$      integer, intent(in) :: iscl
!!$      integer, intent(in), dimension(2) :: bnds
!!$      integer, intent(in) :: vals
!!$      integer :: ix,iy,iz,xi
!!$
!!$      xi=dx*real(bnds(2))
!!$      do iy=iys,iye
!!$         do iz=izs,ize
!!$            do ix=bnds(1),bnds(2)
!!$               t(ix,iy,iscl,iz) = (vals/xi)*(xi-dx*real(ix))
!!$            end do
!!$         end do
!!$      end do
!!$
!!$    end subroutine

!!$    subroutine hygradsource(iscl, bnds, vals)
!!$      integer, intent(in) :: iscl
!!$      integer, intent(in), dimension(2) :: bnds
!!$      integer, intent(in) :: vals
!!$      integer :: ix,iy,iz,yi
!!$
!!$      yi=dy*real(bnds(2))
!!$      do iy=bnds(1),bnds(2)
!!$         do iz=izs,ize
!!$            do ix=1,nnx
!!$               if ((iy >= iys) .and. (iy <= iye)) then
!!$                  t(ix,iy,iscl,iz) = (vals/yi)*(yi-dy*real(iy))
!!$               endif
!!$            end do
!!$         end do
!!$      end do
!!$
!!$    end subroutine
!!$
!!$    subroutine pointsource(iscl, points, vals)
!!$      integer, intent(in) :: iscl
!!$      integer, intent(in), dimension(3) :: points
!!$      integer, intent(in) :: vals
!!$      integer :: ix,iy,iz
!!$      do iy=points(2),points(2) + 2
!!$         do iz=points(3),points(3) + 2
!!$            do ix=points(1),points(1) + 2
!!$               if ((iz >= izs).and.(iz <= ize) .and. (iy >= iys).and.(iy <= iye)) then
!!$                  t(ix,iy,iscl,iz) = vals
!!$               endif
!!$            end do
!!$         end do
!!$      end do
!!$
!!$    end subroutine

    function znptobnd(zt,np)
      integer, intent(in) :: zt
      integer, intent(in) :: np
      integer, dimension(2) :: znptobnd
      integer :: iz

      ! set the first bound, and make sure it doesn't exceed dimensions
      iz = ztoiz(zt)
      znptobnd(1) = iz - int((np-1)/2)
      if (znptobnd(1) < 0) then
        znptobnd(1) = 0
      end if

      ! set the second bound based upon the first
      znptobnd(2) = znptobnd(1) + np -1

    end function

    function xnptobnd(xt,np)
      integer, intent(in) :: xt
      integer, intent(in) :: np
      integer, dimension(2) :: xnptobnd
      integer :: ix

      ! set the first bound, and make sure it doesn't exceed dimensions   
      ix = xtoix(xt)
      xnptobnd(1) = ix - int((np-1)/2)
      if (xnptobnd(1) < 0) then
        xnptobnd(1) = 0
      end if

      ! set the second bound based upon the first                         
      xnptobnd(2) = xnptobnd(1) + np -1

    end function

    function ynptobnd(yt,np)
      integer, intent(in) :: yt
      integer, intent(in) :: np
      integer, dimension(2) :: ynptobnd
      integer :: iy

      ! set the first bound, and make sure it doesn't exceed dimensions   
      iy = ytoiy(yt)
      ynptobnd(1) = iy - int((np-1)/2)
      if (ynptobnd(1) < 0) then
        ynptobnd(1) = 0
      end if

      ! set the second bound based upon the first                         
      ynptobnd(2) = ynptobnd(1) + np -1

    end function

    function restobnd(zt,dr)
      integer,intent(in) :: zt
      integer, intent(in) :: dr
      integer, dimension(2) :: restobnd
      integer :: iz

      iz = ztoiz(zt)
      if (dr > 0) then ! surface res
        restobnd(1) = 0
        restobnd(2) = iz
      else
        restobnd(1) = 0
        restobnd(2) = nnz
      end if
    end function

    function ztoiz(zt)
      integer, intent(in) :: zt
      integer :: ztoiz

      ! note that this will only work for equispaced z grids
      ztoiz = int(zt/dz)

    end function
    
    function xtoix(xt)
      integer, intent(in) :: xt
      integer :: xtoix

      ! note that this will only work for equispaced z grids              
      xtoix = int(xt/dx)

    end function

    function ytoiy(yt)
      integer, intent(in) :: yt
      integer :: ytoiy

      ! note that this will only work for equispaced z grids              
      ytoiy = int(yt/dy)

    end function

end module

