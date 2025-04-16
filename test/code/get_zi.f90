SUBROUTINE get_zi(gradmax,gradout,len,itype)

  USE pars

  REAL :: gradmax(*), gradout(*)

  DO i=1,len,2
    IF(gradmax(i) > gradout(i)) THEN
      gradout(i)   = gradmax(i)
      gradout(i+1) = gradmax(i+1)
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE
