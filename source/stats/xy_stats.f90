SUBROUTINE midpnt(a,b,s,n)

  INTEGER it,j

  IF(n .eq. 1) THEN
    s = (b - a)*stokes_ker(0.5*(a+b))
  ELSE
    it   = 3**(n-2)

    tnm  = FLOAT(it)
    del  = (b - a)/(3.0*tnm)
    ddel = del + del
    x    = a + 0.5*del
    sum = 0.0

    DO j=1,it
      sum = sum + stokes_ker(x)
      x   = x + ddel
      sum = sum + stokes_ker(x)
      x   = x + del
    ENDDO
    s = (s + (b - a)*sum/tnm)/3.0
  ENDIF

  RETURN

! FORMAT
6000  FORMAT(' n = ',i4,' it = ',i4)
1200  FORMAT(' 1200 a = ',e15.6,' b = ',e15.6,/,'      s = ',e15.6,' n = ',i5)

END
