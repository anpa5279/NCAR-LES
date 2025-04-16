SUBROUTINE busngr(zeta,phim,phis,psim,psis)
!BUSINGER VERSION OF SIMILARITY THEORY

  DATA pih /1.57079633/ !pi/2
  SAVE pih

  IF(zeta < 0.) THEN
    x=(1.0 - 15.0*zeta)**0.25
    phim = 1.0/x
    psim = 2.0*ALOG((1.0+x)/2.0) + ALOG((1.0+x*x)/2.0) - 2.0*ATAN(x)+pih

    IF(psim>2.0)psim=2.0

    y = SQRT(1.0-9.0*zeta)
    phis = 0.74/y
    psis = ALOG((1.0+y)/2.0)*2.0
  ELSE IF(zeta > 0) THEN 
    phim = 1.0 + 4.7*zeta
    phis = 0.74 + 4.7*zeta
    psim = -4.7*zeta
    psis = -4.7*zeta
  ELSE
    phim = 1.0
    phis = 0.74
    psim = 0.0
    psis = 0.0
  ENDIF

  RETURN
END SUBROUTINE
