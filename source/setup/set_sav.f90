SUBROUTINE set_sav(it, istart)

    USE pars
    USE fields
    USE con_DATA
    USE con_stats

    DATA ionce/0/
    SAVE ionce

    IF (it /= istart) THEN
        ! INCREMENT TIME IF NOT FIRST TIME THROUGH
        time = time + dt
    END IF

    it = it + 1
    it_counter = it - iti

    dt = dt_new
    dtg = dt
    mnout = (MOD(it_counter, imean) == 0) .OR. (it == 1)
    mtape = (MOD(it_counter, itape) == 0) .OR. (it == 1)
    micut = (MOD(it_counter, itcut) == 0)
    mviz = (MOD(it_counter, i_viz) == 0)

    IF (ihst < 0) THEN
        mhis = .FALSE.
    ELSE
        mhis = (MOD(it_counter, ihst) == 0 .AND. it >= it_his)
    END IF

    ! DECIDE WHETHER VELOCITY FIELDS ARE SAVED
    msave = .FALSE.
    IF (it_counter >= itstr .AND. mtape) THEN
        itn = itn + 1
        msave = .TRUE.
        CALL get_output_filenames
    END IF

    ! DECIDE WHETHER VIZ FIELDS ARE SAVED
    msave_v = .FALSE.
    IF (it_counter >= itstr .AND. mviz .AND. i_viz > 0) THEN
        msave_v = .TRUE.
        IF (ionce == 0) THEN
            ionce = 1
            CALL open_viz
        END IF
    END IF

    ! DECIDE WHETHER HISTORY FILES ARE TO BE SAVED
    IF ((ihst > 0) .AND. (it_counter >= it_nxt)) THEN
        CALL open_his(it)
        it_nxt = it_nxt + itape
    END IF

    RETURN
END SUBROUTINE
