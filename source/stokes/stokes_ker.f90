FUNCTION stokes_ker(sigma)
! EVALUATE KERNEL OF THE STOKES INTEGRAL FOR THE DONELAN SPECTRAL SHAPE
! THE LATTER IS CONVERTED INTO 'SIGMA' SPACE
! CAREFUL WITH 2PI FACTORS

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  wave_spec =  (ann*gravity*gravity/(sigma_p*sigma**4))*EXP(-bnn*(sigma_p/sigma)**4)
  stokes_ker = 2.0*(wave_spec*sigma**3)*EXP(2.0*sigma*sigma*z_pt/gravity)/gravity

  RETURN
END
