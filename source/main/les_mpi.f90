 PROGRAM les_mpi

  USE pars
  USE fields
  USE con_data
  USE con_stats
  USE tracerbc, ONLY: applytracerbc

  INCLUDE 'mpif.h'

  !INITIALIZE MPI, GET MYID, NUMPROCS
  !TEST IF ON ROOT PROCESS


  !Create child processes, each of which has its own variables. From this point on, every process executes a separate copy of this program.  
  !Each process has a different process ID, ranging from 0 to num_procs minus 1, and COPIES of all variables defined in the program. 
  !No variables are shared.
  CALL mpi_init(ierr) !causes creation of additional processes on different systems and continue executing the program. initializes MPI enviroment

  !find out MY process ID, and how many processes were started.
  CALL mpi_comm_rank(mpi_comm_world,myid,ierr) !myid= which process it is OR process rank 
  CALL mpi_comm_size(mpi_comm_world,numprocs,ierr) !numprocs= total number of processes OR size of cluster

  i_root = 0
  l_root = .FALSE.
  IF(myid == i_root) l_root = .TRUE. !if it is the first process

  l_debug = .FALSE.
  IF(idebug == 1) l_debug = .TRUE. !if debug is on debugging is true

  ts_mpi = mpi_wtime() !mpi_wtime returns the time in seconds since an some time stamp in the past

  !SET NUMBER OF X-Y SLAB CPUS (aka, . this directly relates to the parallelization of the cpus)
  ncpu_s   = 8 !
  itn      = 0 !the counter for the time step saving
  case_inp = '30L'

  CALL get_units !used to assign unit numbers to different files and open a unit for standard print output.
  CALL gridd!is responsible for allocating memory for various arrays used in a computational fluid dynamics (CFD) simulation. 
  CALL setcon
  CALL set_paths
  !istop = 1 !this only shows up here

  !SCRATCH RUN
  IF (iti==0)  THEN !iti defined in pars.f90
    igrdr = 2 !data comes from initialization (random)
    case = case_inp
    CALL init !only time this is ever called/used
    CALL setup(it) !this is where the initial magnitude is defined

    !CHOOSE ROUTINE FOR GETTING INITIAL GUESS
    CALL randoc !only time this is ever called/used
    CALL get_max
    CALL applytracerbc(it) !defining initial conditions for scalars
  ELSE
    igrdr = 3 !data comes from restart file
    CALL restart
    CALL get_max
    CALL setup(it)
  ENDIF

  !TIME LOOP
  tzero = time
  CALL get_dt(it,iti) !setup based on CFL number 0.50 (either fixed time step or adaptive. CFL set in setup/estcon)

  DO WHILE(it<itmax) !itmax is set in pars.f90, it= iterations (dt is the magnitude of timestep)
    CALL set_sav(it,iti) !updates iteration it=it+1 and decides whether to save as output or not

    !UPDATE POSITION OF VORTEX (this if statement seems to be true for every time step unless ivis0 does not equal 1)
    IF(it >= new_vis .AND. ivis0 == 1) THEN 
      ivis = 1 !updated vortex position
    ELSE
      ivis = 0 !vortex position is not updated
    ENDIF

    !3 STAGE RUNGE-KUTTA TIME STEPPING
    DO istage=1,3 !DO 8999 istage=1,3 !looking at RK3 manipulation terms to update time step
      dtzeta = dt*zetas(istage) !units= time. used in RK3 solving solve/comp1
      dtgama = dt*gama(istage) !units= time. used in comp_p and comp2

      !COMPUTE DERIVATIVES OF (U,V,W)
      CALL exchange !mpi process
      CALL get_derv !gets derivatives with respect to x and y

      !NEW EDDY VISCOSITY, AND BCS
      IF(iss == 0 .AND. ifree == 0) THEN !iss=starting processor, ifree=flag (when ifree=0, use spatially averaged surface conditions for MO )
        CALL lower(it) !not communicating between levels. will not be called because it is the lower boundary 
      ELSEIF(ifree == 1) THEN !use point-by-point conditions for MO free convection 
        CALL lower_free(it)
      ENDIF

      IF(ise == numprocs-1) THEN !ise= ending processor, if ending processor equals the size of the cluster minus 1, then last boundary to work on
        CALL upper !you are therefore looking at upper boundary 
      ENDIF

      CALL bcast_pbc !SEND UPPER BC TO OTHER PROCESSORS FOR FFT SOLUTION OF PRESSURE
      CALL get_means(istage) !GET MEANS FOR ALL VARIABLES FOR USE IN ISO, SURFVIS, COMP1, COMPMN

      IF(ivis == 1) THEN !if the vortex position has been updated
        CALL iso(it) !look into isotropy
        CALL surfvis(it)
      ENDIF

      IF(istage == 1)THEN !looking at 1st RK3 manipulation term
        CALL xy_stats !GET STATISTICS
        CALL tke_budget !GET TERMS IN RESOLVED SCALE TKE BUDGET AS IN GABLS WRITEUP AT W-POINTS AT
        CALL pbltop(itop) !getting estimate for planetary boundary layer according to method type (methods described in misc/pbltop)
      ENDIF

      !GET RHS FOR ALL EQUATIONS
      IF(istage==1 .AND. flg_reaction==1)THEN !looking at 1st RK3 manipulation term and if reactions occur (before advection step). the code can be faster if this was nested in the previous if
        CALL strang1(it) !STRANG SPLITTING OF SCALAR REACTION TERM
      ENDIF

      CALL comp1(istage,it) !3RD ORDER RK (RK3) TIME STEPPING AND MONOTONE SCALAR FLUXES IN X,Y,Z DESIGNED TO USE 

      !SOLVE FOR PRESSURE
      CALL comp_p !SETUP PRESSURE SOLVER only 

      !ADD PRESSURE GRADIENT AND DEALIAS
      CALL comp2 !solving for rhs and updating stuff for RK3

      IF(istage==3 .AND. flg_reaction==1)THEN !looking at 3rd RK3 manipulation term and if reactions occur (after the advection step)
        CALL strang1(it) !STRANG SPLITTING OF SCALAR REACTION TERM
      ENDIF


      !saving variables
      IF(msave .AND. istage == 3) THEN !if final step of RK3 and if msave is TRUE
        CALL save_v(it) !save 3D field
      ENDIF

      IF(istage == 3) THEN !we can clean this up by doing nested if statements
        IF(msave .AND. l_root) CALL save_c(it)
      ENDIF


      IF(micut) THEN !micut = (MOD(it_counter,itcut)==0)
        CALL dealias !wave cutoff filter using 2d fft
      ENDIF

      IF(mnout .AND. istage == 1)  THEN !mnout = (MOD(it_counter,imean)==0), no remainder .OR. (it==1), first time step that is not the initial condition
        IF(l_debug) THEN !if debugging print
          CALL print(nprt,it,izs,ize)
        ENDIF
        IF(l_root) CALL print(6,it,1,nnz)
      ENDIF

      IF(l_root) THEN
        IF(mhis  .AND. istage == 1)  CALL write_his(itop)
        IF(mhis  .AND. istage == 1 .AND. mtape) CALL close_his
      ENDIF

    END DO !8999 CONTINUE !marks the end of the do loop (apparently outdated)
    CALL get_max
    CALL get_dt(it,iti)
    
  END DO !completes the do while loop


  te_mpi = mpi_wtime() !mpi_wtime returns the time in seconds since an some time stamp in the past (since line 31)

  WRITE(6,9997) (te_mpi - ts_mpi) !write the time difference from start to finish in seconds

  CONTINUE
  CALL mpi_finalize(ierr) !completes mpi_init and ends MPI communications

!FORMAT
9997  FORMAT(' Job Execution Time = ',e15.6)

  STOP
END PROGRAM
