SUBROUTINE get_units
!used to assign unit numbers to different files and open a unit for standard print output.
  USE pars

  ! UNIT NUMBERS FOR FILES
  nvel  = 20
  npre  = 30
  nhis1 = 40
  nvelc = 50 !this is the number for the files in the data folder
  nhisp = 60
  nviz_z = 80
  nviz_y = 82
  nviz_x = 84
  nviz_s = 90

  ! UNIT NUMBER FOR STANDARD PRINT OUT FOR EACH MPI TASK
  nprt = 1

  ! OPEN UNIT FOR STANDARD PRINT OUT
  path_prt = case_inp(1:3)//'.mp.xxxxx.out' !creating the files that are outside the data folder
  WRITE(path_prt(8:12),'(i5.5)') myid !myid= which process this is
  OPEN(nprt,file=path_prt,form='formatted')

  RETURN
END SUBROUTINE
