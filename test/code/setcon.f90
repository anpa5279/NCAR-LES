SUBROUTINE setcon

  USE pars
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'
  EXTERNAL get_zi

  !GET MACHINE TYPE, CAN ERAD DATA DIRECTORY
  !INITIALIZE FFT
  CALL rffti(nnx,trigx(1,1))
  CALL rffti(nny,trigx(1,2))
  CALL cffti(nny,trigc(1))

  !START STEP FOR HISTORY FILES
  it_nxt = it_his

  !SET MIN VALUE OF E
  smal_e = 0.0
  smal_e = 1.0e-12

  !SET CONSTANTS IN EDDY VISCOSITY MODEL
  ck       = 0.1 !open channel flow
  ceps     = 0.93
  csmag    = SQRT(ck*SQRT(ck/ceps))
  stab_c   = 0.76

  !SET STABILITY CONSTANT
  stabmin = 1.0e-12

  !MINIMUM DSL LENGTH CONSTANT
  almin_c = 0.0001

  !INITIALIZE GRID RESTART FLAG
  igrdr = 1

  !CREATE MPI OPERATION TO FIND MAX, AND LOCATION USING LOCAL GRADIENT METHOD
  CALL mpi_op_create(get_zi,.TRUE.,ziloc,ierror)

  !DEFINE COEFFS FOR O(3) RK TIME STEPPING SCHEME, BURROWED FROM SPALART,
  !MOSER, AND ROGERS (JCP, 1990)
  !THIS IS A SIMPLER VERSION SINCE ALL TERMS ARE LUMPED IN NON-LINEAR TERMS.
  !CFL NUMBER IS FOR AN ENTIRE RK STEP IN THE CASE THREE STAGES
  !CFL = MAX(U)*DT/DX
  !these are rates
  zetas(1) = 0.0
  zetas(2) = -17.0/60.0
  zetas(3) = -5.0/12.0
  gama(1)  = 8.0/15.0
  gama(2)  = 5.0/12.0
  gama(3)  = 3.0/4.0
  !etas(1)  = -1.0
  !etas(2)  = -1.0 + 8.0/15.0
  !etas(3)  = -1.0 + 2.0/3.0

  !A FULL STEP, AT THE NEW TIME
  !etas(4)  =  0.0

  cfl = 0.50

  RETURN
END SUBROUTINE
