SUBROUTINE lterp(n,zary,zpt,i,ip1,ratio)
!LINEAR INTERPOLATION FOR ZPT IN ZARY, WHERE ZARY IS MONOTONIC INCREASING OR
!DECREASING FUNCTION

  DIMENSION zary(*)
  nm1 = n-1

  IF(n<=1) THEN
    i = 1
    ip1 = 1
    ratio = 0.0
    go to 999
  ENDIF

  IF(zary(1) < zary(2)) go to 1
     go to 101
  1 CONTINUE

  !MONOTONIC INCREASING ARRAY
  IF(zpt < zary(1)) THEN
    i = 1
    ip1 = 1
    ratio = 0.0
    go to 999
  ELSE IF(zpt > zary(n)) THEN
    i = n
    ip1 = n
    ratio = 1.0
    go to 999
  ENDIF

  DO j=1,nm1
    IF(zpt >= zary(j) .AND. zpt <= zary(j+1)) THEN
      i = j
      ip1 = j+1
      ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
      go to 999
    ENDIF
  ENDDO

  !DECREASING ARRAY
  101 CONTINUE

  IF(zpt > zary(1)) THEN
    i = 1
    ip1 = 1
    ratio = 0.0
    go to 999
  ELSE IF(zpt < zary(n)) THEN
    i = n
    ip1 = n
    ratio = 1.0
    go to 999
  ENDIF

  DO j=1,nm1
    IF(zpt <= zary(j) .AND. zpt >= zary(j+1)) THEN
      i = j
      ip1 = j+1
      ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
      go to 999
    ENDIF
  ENDDO

  999 CONTINUE
  RETURN
END SUBROUTINE
