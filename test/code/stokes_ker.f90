FUNCTION stokes_ker(sigma)
!EVALUATE KERNEL OF THE STOKES INTEGRAL FOR THE DONELAN SPECTRAL SHAPE
!THE LATTER IS CONVERTED INTO 'SIGMA' SPACE
!CAREFUL WITH 2PI FACTORS

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  wave_spec =  (ann*grav*grav/(sigma_p*sigma**4))*EXP(-bnn*(sigma_p/sigma)**4) !eq 33 in Webb and Fox-Kemper 2011 without gamma_DHH and Gamma_DVH
  stokes_ker = 2.0*(wave_spec*sigma**3)*EXP(2.0*sigma*sigma*z_pt/grav)/grav !eq 41 in Webb and Fox-Kemper 2011 and eq 43 but the integral is not computed.
  !z_pt is defined in stokesv

  RETURN
END
