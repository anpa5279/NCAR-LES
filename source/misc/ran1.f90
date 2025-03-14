FUNCTION ran1(idum)
! STOLEN FROM NUMERICAL RECIPES, P. 271

  INTEGER, PARAMETER :: ntab = 32
  INTEGER :: idum, ia, im, iq, ir, ndiv
  REAL :: ran1, am, eps, rnmx
  INTEGER :: j, k, iv(ntab), iy

  SAVE iv, iy
  DATA iv /ntab*0/, iy /0/

  ia=16807
  im=2147483647
  iq=127773

  am=1.0/im
  ir=2836.0
  ndiv=1+(im-1)/ntab
  eps=1.2e-07
  rnmx=1.0-eps

  IF(idum <= 0 .OR. iy == 0) THEN
    idum = MAX(-idum,1)
    DO j=ntab+8,1,-1
      k = idum/iq
      idum = ia*(idum - k*iq) - ir*k
      IF(idum < 0) idum = idum + im
      IF(j <= ntab) iv(j) = idum
    ENDDO
    iy = iv(1)
  ENDIF

  k     = idum/iq
  idum  = ia*(idum - k*iq) - ir*k

  IF(idum < 0) idum = idum + im

  j     = 1 + iy/ndiv
  iy    = iv(j)
  iv(j) = idum
  ran1  = MIN(am*iy, rnmx)

  RETURN
END
