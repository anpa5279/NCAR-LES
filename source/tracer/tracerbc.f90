module tracerbc
  use con_data, only: dz,dx,dy
  use con_stats, only: z, zz
  use pars, only: iys,iye,izs,ize,izi,nnx,nnz,nny,nscl
  use fields, only: t
  use inputs

  implicit none
  integer, parameter :: flg_debug = 0
  real, dimension(nscl) :: tau,airval
  integer, dimension(nscl) :: ictype,rmodel,rdorg,rpartner,asflux
  integer, dimension(2,nscl) :: bnd
  integer, dimension(2) :: bnds
  integer, dimension(3,nscl) :: point
  integer, dimension(3) :: points
  real, dimension(nscl) :: val
  contains

! iscl      : scalar number (temperature is always iscl=1)
! tau       : reaction time scale
! ictype    : initial condition (0 = nothing, 1 = horiz. band,
!                                2 = vertical band in x, 3 = vertical band in y,
!                                4 = point source, 5 = vertical gradient,
!                                6 = horiz. gradient in x, 7 = horiz. gradient in y)
!             ictype does not work for iscl=1 (temperature), that is set in init/randoc.f
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
! asflux    : air-sea flux boundary condition (0 = for no flux, 1 = for flux [also need
!             flag_airseaflux.eq.1 in pars.f])
! airval    : value of tracer in air (only need set for use with asflux and flag_airseaflux)

    subroutine applytracerbc(it)
      integer, intent(in) :: it
      integer :: iscl, np, zt
      real :: ta, vals

      !! active tracers (temperature)
            iscl = 1;
            ictype(iscl) = 0;   val(iscl) = 273.15 + iTsurf;      tau(iscl)    = 0;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = 0;      zt = 0;  rmodel(iscl) = 0;  bnd(:,iscl) = znptobnd(zt,np);

            !! passive tracers
            iscl = 2; !carbonate dioxide
            ictype(iscl) = 1;   val(iscl) = c1;     tau(iscl)      = 1;
            asflux(iscl) = 1;   airval(iscl) = 8.56056;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 3; !bicarboante
            ictype(iscl) = 1;   val(iscl) = c2;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 4; !carbonate
            ictype(iscl) = 1;   val(iscl) = c3;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 5; !boric acid
            ictype(iscl) = 1;   val(iscl) = c4;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 6; !tetrahydroborate
            ictype(iscl) = 1;   val(iscl) = c5;  tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 7; !hydrogen ion
            ictype(iscl) = 1;   val(iscl) = c6; tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            iscl = 8; !hydroxyl
            ictype(iscl) = 1;   val(iscl) = c7;     tau(iscl)      = 1;
            asflux(iscl) = 0;   airval(iscl) = 0;
            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);
!=================NPZD========================
            !            iscl = 9; !p
            !            ictype(iscl) = 1;   val(iscl) = c8;  tau(iscl)      = 1;
            !            asflux(iscl) = 0;   airval(iscl) = 0;
            !            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            !            iscl = 10; !z
            !            ictype(iscl) = 1;   val(iscl) = c9; tau(iscl)      = 1;
            !            asflux(iscl) = 0;   airval(iscl) = 0;
            !            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

            !            iscl = 11; !n
            !            ictype(iscl) = 1;   val(iscl) = c10;     tau(iscl)      = 1;
            !            asflux(iscl) = 0;   airval(iscl) = 0;
            !            np = nnz+2;  zt = 0;  rmodel(iscl) = 3;  bnd(:,iscl) = znptobnd(zt,np);

      do iscl = 2,nscl !skipping temperature
        bnds=bnd(:,iscl); vals=val(iscl); points=point(:,iscl);
        if (ictype(iscl).eq.1) call hbndsource(iscl,bnds,vals); 
        if (ictype(iscl).eq.5) call vgradsource(iscl,bnds,vals);
        bnds = 0; vals = 0; points = 0;
      enddo

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

    subroutine vgradsource(iscl, bnds, vals)
      integer, intent(in) :: iscl
      integer, intent(in), dimension(2) :: bnds
      real, intent(in) :: vals
      integer :: ix,iy,iz,zi

      zi  = z(bnds(2))
      do iy=iys,iye
         do iz=bnds(1),bnds(2)
            do ix=1,nnx
               if ((iz >= izs) .and. (iz <= ize)) then
                  t(ix,iy,iscl,iz) = (vals/zi)*(zi-zz(iz))
               endif
            end do
         end do
      end do

    end subroutine

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
