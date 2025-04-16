SUBROUTINE set_sav(it,istart)

  USE pars
  USE fields
  USE con_DATA
  USE con_stats

  DATA ionce /0/
  SAVE ionce

  IF(it /= istart) THEN
    !INCREMENT TIME IF NOT FIRST TIME THROUGH
    time=time+dt
  ENDIF

  it=it+1
  it_counter = it - iti

  dt    = dt_new !updates new time step 
  dtg   = dt
  mnout = (MOD(itcounter,imean)==0) .OR. (it==1)
  mtape = (MOD(it_counter,itape)==0)
  micut = (MOD(it_counter,itcut)==0)
  mviz  = (MOD(it_counter,i_viz)==0)

  IF(ihst < 0) THEN
    mhis = .FALSE.
  ELSE
    mhis = (MOD(it_counter,ihst)==0 .AND. it >= it_his)
  ENDIF

  !DECIDE WHETHER VELOCITY FIELDS ARE SAVED
  msave = .FALSE.
  IF(it_counter >= itstr .AND. mtape) THEN !if it_counter>1 and it can be divided by itape without a remainder
    itn=itn+1
    msave = .TRUE.
    CALL get_output_filenames
  ENDIF

  !DECIDE WHETHER VIZ FIELDS ARE SAVED
  msave_v = .FALSE.
  IF(it_counter >= itstr .AND. mviz .AND. i_viz > 0) THEN
    msave_v = .TRUE.
    IF(ionce == 0) THEN
      ionce = 1
      CALL open_viz
    ENDIF
  ENDIF

  !DECIDE WHETHER HISTORY FILES ARE TO BE SAVED
  IF((ihst > 0) .AND. (it_counter >= it_nxt)) THEN
    CALL open_his(it)
    it_nxt = it_nxt + itape
  ENDIF

  RETURN
END SUBROUTINE
