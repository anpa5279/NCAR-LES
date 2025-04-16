SUBROUTINE strang1(it)
!STRANG SPLITTING OF SCALAR REACTION TERM - 0.5*REACT, ADVECT, 0.5*REACT FOR
!FAST REACTIONS (TAU <= 1000), SEE RHS/RHS_SCL.F90 FOR SLOW REACTION SOURCES
!Strang splitting is a numerical method for solving differential equations that are decomposable into a sum of differential operators. used to speed up calculation for problems involving operators on very different time scales, for example, chemical reactions in fluid dynamics, and to solve multidimensional partial differential equations by reducing them to a sum of one-dimensional problems
!only called in les_mpi when looking at 1st and 3rd RK3 manipulation term and if reactions occur
!chemistry is integrated in two half steps, before and after the advection step
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats
  USE reaction, ONLY: react_src !only time reaction is used

  INCLUDE 'mpif.h'

  !TEMP SCALAR ARRAY TO HOLD RHS FROM STEP N-1 AND FIELD VARIABLES FROM STEP N
  REAL :: trhs(nnx,iys:iye,nscl,izs:ize)
  REAL, DIMENSION(nscl-1) :: tmp

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
      tmp = react_src(ix,iy,1,iz) !that 1 in the code does not do anything
        DO l=2,nscl
          trhs(ix,iy,l,iz) = tmp(l-1)
          IF(trhs(ix,iy,l,iz).le.1.0e-20)THEN
            trhs(ix,iy,l,iz) = 1.0e-20
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO iz=izs,ize
    DO l=2,nscl
      DO iy=iys,iye
        DO ix=1,nnx
          t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
