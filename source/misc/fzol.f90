SUBROUTINE fzol(zeta,phim,phis,psim,psis)
! ESTIMATE THE STABILITY FUNCTION FOR MOMENTUM, M AND SCALARS, C FROM INPUT OF
! THE STABILITY PARAMETER ZETA = Z/L

  DATA c1/5./
  DATA a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
  DATA zetam,zetas/-0.2,-1.0/
  SAVE c1, a3, b3, a4, b4, zetam, zetas

  psimu(Y)  = 1.571 + 2.0*(ALOG(0.5*(1.0 + Y)) - ATAN(Y)) + ALOG(0.5 + 0.5*Y**2)
  psisu(Y)  = 2.0*ALOG(0.5 + 0.5*Y)
  psicu(Y,G)= (1.0 - G)*ALOG(ABS(Y - 1.0)) + 0.5*(G + 2.0)*ALOG(ABS(Y**2 +  &
        Y + 1.0)) - (2.0*G + 1.0) / SQRT(3.0) * ATAN((Y + 0.5)*2.0/SQRT(3.0))

  Xm(zol)   = (1.0 - 16.0*zol)**0.25
  Xs(zol)   = SQRT(1.0 - 16.0*zol)
  Xc(zol,f) =  ABS(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)

  IF(zeta>=0.0) THEN
    IF(zeta<=1.0) THEN
      phim = 1.0 + c1 * zeta
      psim = - c1 * zeta
      phis = phim
      psis = psim
    ELSE
      ! USE LIMITING FORM
      phim = c1 + zeta
      psim = (1.0 - c1)*(1.0 + ALOG(zeta) ) - zeta
      phis = phim
      psis = psim
    ENDIF
  ELSE
    IF(zeta.ge.zetam) THEN
      phim = 1.0 / Xm(zeta)
      psim = psimu(Xm(zeta))
    ELSE
      ! USE CONVECTIVE LIMIT FOR MOMENTUM
      X = (1.0 - b3/a3 * zeta)**(1.0/3.0)
      fm = a3**(-1.0/3.0)
      phim = fm / Xc(zeta,b3/a3)
      psim = psimu(Xm(zetam)) + psicu(Xc(zeta,b3/a3),fm) -                  &
            psicu(Xc(zetam,b3/a3),fm)
    ENDIF

    ! UNSTABLE SCALARS
    IF(zeta.ge.zetas) THEN
      phis = 1.0/Xs(zeta)
      psis = psisu(Xs(zeta))
    ELSE
      ! USE CONVECTIVE LIMIT FOR SCALARS
      fs =   ABS(a4)**(-1.0/3.0)*ABS(a4)/a4
      phis = (a4 - b4*zeta)**(-1.0/3.0)
      psis = psisu(Xs(zetas))+psicu(Xc(zeta,b4/a4),fs) - psicu(Xc(zetas,b4/a4),fs)
    ENDIF

  ENDIF
  RETURN
END SUBROUTINE
