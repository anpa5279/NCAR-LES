SUBROUTINE stokesv
! GET STOKES DRIFT VELOCITY USING DONELAN SPECTRUM MATCHED TO WAVE DATA
! SEE ALVES ET AL., JPO 2003, VOL. 33
! COMPARE ALL PDF VARIABLES IN SR. PDF_NDOT
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'

  ! DONELAN SPACE (CONSTANTS SET IN INIT)
  ! RECOMPUTE PEAK IF WE CHANGE U_10 !!!
  CALL speed2stress(u_10,v_10,cd_10,tau_x,tau_y)

  speedval = SQRT(u_10**2 + v_10**2)
  f_p     = f2w*grav/speedval
  sigma_p = two_pi*f_p

  ! SET PARAMETERS THOUGH NOT USED HERE
  stokesw = 0.0
  stokesa = 0.0
  stokess = 0.0

  range_min = 0.1
  range_max = 5000.0
  DO iz=1,nnzp1
    z_pt = zz(iz)
    CALL s_int(range_min,range_max,value)
    stokes(iz) = value

    ! FOR NO STOKES
    IF(flg_stokes /= 1) THEN
      stokes(iz) = 0.0
    ENDIF
  ENDDO

  IF(l_root) THEN
    WRITE(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
  ENDIF

  dir_x = 1.0
  dir_y = 0.0

  RETURN

! FORMAT
2212  FORMAT('#k ',/,'#m 4',/,'#lw 1.0',/,(2e15.6))
6000  FORMAT(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))

END SUBROUTINE
