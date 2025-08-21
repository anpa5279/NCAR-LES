SUBROUTINE speed2stress(u_10,v_10,cd_10,tau_x,tau_y)

  USE pars
  USE inputs

  ! SPECIFY SPEED HERE (KEEP BETWEEN 5-10 M/S) -- why???
  IF(flg_lat /= 1) THEN
    u_10 = ustokes
  ELSE
    !u_10 = SQRT(taux_est/(cd_fac*(DTANH((SQRT(taux_est/(1.3e-3)) - 25.0)/5.0)*0.5 - DTANH(-25.0/5.0)*0.5 + 1.2)*0.001))
  ENDIF

    v_10 = 0.0
    s_10 = SQRT(u_10**2 + v_10**2)

    arg1   = (s_10 - 25.0)/5.0
    arg2   = -25.0/5.0
    a1     = TANH(arg1)
    a2     = TANH(arg2)
    cd_10  = cd_fac*(a1*0.5 - a2*0.5 + 1.2)*0.001

    ! LIMITED CD
    tau_x  = cd_10*s_10*u_10 !this does not manipualte anything 
    tau_y  = cd_10*s_10*v_10

  RETURN
END SUBROUTINE


! could i simplify the coriolis effects out for the purposes of this study and
! simplifying the stokes drift calculation?
