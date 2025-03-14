SUBROUTINE s_int(r_min,r_max,value)
!GET INTEGRAL USING A MID-POINT RULE STOLEN FROM NUMERICAL RECIPES

  iter = 10
  value = 0

  DO j=1,iter
    CALL midpnt(r_min,r_max,value,j)
  ENDDO

  RETURN

!FORMAT
1000  FORMAT(' j = ',i5,' value = ',e15.6)

END SUBROUTINE
