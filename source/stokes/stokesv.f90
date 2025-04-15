SUBROUTINE stokesv
!GET STOKES DRIFT VELOCITY USING DONELAN SPECTRUM MATCHED TO WAVE DATA
!SEE ALVES ET AL., JPO 2003, VOL. 33
!COMPARE ALL PDF VARIABLES IN SR. PDF_NDOT
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'

  !DONELAN SPACE (CONSTANTS SET IN INIT) (Donelan et al. 1985)
  !RECOMPUTE PEAK IF WE CHANGE U_10 !!!
  CALL speed2stress(u_10,v_10,cd_10,tau_x,tau_y) !these variables are only used in speed2stress.f90

  speedval = SQRT(u_10**2 + v_10**2) !v_10 is 0, so this is just u_10
  f_p     = f2w * grav / speedval !does not change 
  sigma_p = 2 * pi * f_p ! does not change 

  !SET PARAMETERS THOUGH NOT USED HERE
  stokesw = 0.0
  stokess = 0.0

  range_min = 0.1
  range_max = 5000.0
  DO iz=1,nnzp1
    z_pt = zz(iz) !middle of the cell
    CALL s_int(range_min,range_max,value) !value = varying stokes vertically with depth. this calls midpnt which uses function stokes_ker to solve the vertical direction. 
    stokes(iz) = value
  ENDDO

  IF(l_root) THEN
    WRITE(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
  ENDIF

  dir_x = 1.0
  dir_y = 0.0

  RETURN

!FORMAT
2212  FORMAT('#k ',/,'#m 4',/,'#lw 1.0',/,(2e15.6))
6000  FORMAT(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))

END SUBROUTINE
