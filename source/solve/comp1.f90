SUBROUTINE comp1(istep,it)
  !called in les_mpi always
!3RD ORDER RK (RK3) TIME STEPPING AND MONOTONE SCALAR FLUXES IN X,Y,Z DESIGNED TO USE
!MPI IN X AND Y DIRECTIONS (aka designed to work with the x-y layers that span in the z-direction)
!istep = istage = which step of RK3
!it = iterations (the time stepping aspect)
!T_DISS AND TR_TAU ARE IN COMP1
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h' !this implements mpi for the Fortran interface. All program units that make MPI calls must include this. apparently discouraged and may be gone in future MPI versions

  INTEGER :: istatus(mpi_status_size) 
  INTEGER, PARAMETER ::                   &
      js = 7,                            & !NUMBER OF NON-SCALAR STATS 
      ns = 3                                !NUMBER OF SCALAR STATS
  INTEGER, PARAMETER :: nstat = js + ns*nscl
  REAL :: stat(1:nnz,nstat)

  !TEMP ARRAYS TO HOLD RHS FROM STEP N-1 AND FIELD VARIABLES FROM STEP N
  REAL :: urhs(nnx,iys:iye,izs:ize),vrhs(nnx,iys:iye,izs:ize),              &
          wrhs(nnx,iys:iye,izs:ize),erhs(nnx,iys:iye,izs:ize),              &
          trhs(nnx,iys:iye,nscl,izs:ize)

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        !applying RK3
        !rs are velocity defined in rhs/rhs_uvw
        !dtzeta is the time*constant defined in les_mpi
        urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
        vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
        wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
        !r5 is kinetic energy defined in diff/tke_vis
        erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  DO iz=izs,ize !looping z
    DO l=1,nscl !looping scalars. nscl= number of scalars and vars (when nscl=1 then you only observe temperature)
      DO iy=iys,iye !looping y
        DO ix=1,nnx !looping x
        !applying RK3
        !r4 is defined in rhs/rhs_scl
          trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz) !t= scalars
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !GET VISCOSITY AND RHS OF (E,U,V,W) EQUATIONS AT NEXT STEP
  CALL tke_vis(istep) !only called here
  CALL rhs_uvw(istep) !only called here. updating values for RK3

  !EVALUATE RHS OF SCALAR EQUATIONS=
  DO l=1,nscl !looping scalars
    CALL rhs_scl(istep,l,it)!only called here. updating values for RK3 for every scalar
  ENDDO

  !GATHER STAT SUMS ON ROOT PROCESSOR USING MPI_REDUCTION OVER ALL PROCESSORS
  IF(istep == 1) THEN !if in the first step of RK3
    DO j=1,nstat
      DO iz=1,nnz
        stat(iz,j) = 0.0
      ENDDO
    ENDDO

    DO iz=izs,ize
      stat(iz,1) = uwsb(iz)
      stat(iz,2) = vwsb(iz)
      stat(iz,3) = wwsb(iz)
      stat(iz,4) = tr_tau(iz)
      stat(iz,5) = triz(iz)
      stat(iz,6) = shrz(iz)
      stat(iz,7) = t_diss(iz)
    ENDDO

    m1 = js
    m2 = js + nscl !nscl= number of scalars and vars (nscl=1 is only physics)
    m3 = js + 2*nscl
    DO l=1,nscl
      DO iz=izs,ize
        stat(iz,m1+l) = utsb(iz,l)
        stat(iz,m2+l) = vtsb(iz,l)
        stat(iz,m3+l) = wtsb(iz,l)
      ENDDO
    ENDDO

    CALL mpi_sum_z(stat(1,1),i_root,myid,nstat*nnz,1)

    DO iz=1,nnz
      uwsb(iz)   = stat(iz,1)
      vwsb(iz)   = stat(iz,2)
      wwsb(iz)   = stat(iz,3)
      tr_tau(iz) = stat(iz,4)
      triz(iz)   = stat(iz,5)
      shrz(iz)   = stat(iz,6)
      t_diss(iz) = stat(iz,7)
    ENDDO

    DO l=1,nscl
      DO iz=1,nnz
        utsb(iz,l) = stat(iz,m1+l)
        vtsb(iz,l) = stat(iz,m2+l)
        wtsb(iz,l) = stat(iz,m3+l)
      ENDDO
    ENDDO

    DO iz=1,nnz
      buyz(iz) = batag*wtsb(iz,1)
    ENDDO
  ENDIF

  !SAVE OLD RHS IN FIELD VARIABLES FOR RK ADVANCEMENT
  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        u(ix,iy,iz) = urhs(ix,iy,iz)
        v(ix,iy,iz) = vrhs(ix,iy,iz)
        w(ix,iy,iz) = wrhs(ix,iy,iz)
        e(ix,iy,iz) = erhs(ix,iy,iz)
      ENDDO
    ENDDO
  ENDDO

  DO iz=izs,ize
    DO l=1,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
