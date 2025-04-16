SUBROUTINE blnk(word)

  CHARACTER :: word*(*)

  nchar = LEN(word)

  DO j=1,nchar
    word(j:j) = ' '
  ENDDO

  RETURN
END SUBROUTINE
