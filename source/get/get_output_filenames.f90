SUBROUTINE get_output_filenames
! BUILD FILE NAMES FOR VELOCITY, PRESSURE, AND CONSTANTS

  USE pars
  INCLUDE 'mpif.h'

  CHARACTER cgrid*10, num*5

  ! BUILD CHARACTER STRINGS FOR FNAME
  cgrid = '.mp.'
  WRITE(num,'(i5.5)') itn
  iblnk = index(path_sav,' ')

  CALL blnk(path_sav_v)
  CALL blnk(path_sav_p)
  CALL blnk(path_sav_c)

  path_sav_v = path_sav(1:iblnk-1)//'/u'//cgrid(1:4)//case(1:3)//num(1:5)
  path_sav_p = path_sav(1:iblnk-1)//'/p'//cgrid(1:4)//case(1:3)//num(1:5)
  path_sav_c = path_sav(1:iblnk-1)//'/u'//cgrid(1:4)//case(1:3)//num(1:5)//'.con'

  RETURN
END SUBROUTINE
