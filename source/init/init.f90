SUBROUTINE init
!environment properties. the context
!only called in les_mpi line 48 when iti=0 (before the while loop. data comes from random)
  USE pars
  USE inputs
  USE fields
  USE con_data
  USE con_stats

  pi     = 4.0*ATAN(1.0)
  d_to_r = 2.0*pi/360.0 !degrees to radians
  grav   = 9.81 !gravity
  bfac   = 1.0 !buoyancy factor

  IF(ibuoy==0) bfac = 0. !no buoyancy in the model

  !CASE SPECIFIC DATA
  rho_a   = 1.0 !density of air
  rho_w   = 1000.0 !density of water
  t00b    = 5000.0 !1/(thermal expansion cofficient)
  cp_a    = 1.0057e03 !pressure coefficient of air
  cp_w    = 4.20e03 !pressure coefficient of water 
  gcp     = grav/cp_w
  batag   = bfac*grav/t00b

  !SPECIFY STOKES DRIFT PARAMETERS
  cpou10  = 0.6
  turb_la = 0.3 !turbulent langmuir number (La_t=sqrt(friction velocity/stokes drift)). Langmuir turbulence occurs when a surface boundary layer is forced by wind in the presence of surface waves. fully developed flow according to Webb and Fox-Kemper 2011. 
  rlat    = 30

  fcor    = 2.0*2.0*pi*SIN(rlat*d_to_r)/(24.0*3600.0)
  fcor_h  = 0.0
  ugcont  = 0.0
  vgcont  = 0.
  wtsfc(1) = hflux !hflux defined in main/inputs line 7
  qstar(1) = wtsfc(1) !hflux

  !OTHER THERMODYNAMIC VARIABLES FOR HURRICANE TIME VARYING HEAT FLUXES
  c_e    = 0.0015         !EXCHANGE COEFFICIENT FOR MOISTURE
  c_h    = 0.0015         !EXCHANGE COEFFICIENT FOR TEMPERATURE
  rlv    = 2.45e06        !LATENT HEAT OF VAPORIZATION J/KG
  p_surf = 1000.0         !SURFACE PRESSURE IN MILLIBARS
  t_surf = 273.15 + iTsurf  !SURFACE TEMPERATURE IN KELVIN
  q_rel  = 0.8            !RELATIVE HUMIDITY IN % AT 10M HEIGHT
  t_10m  = t_surf - 2.5   !10M TEMPERATURE

  !SATURATION VAPOR PRESSURE (GILL, P606). eventually used in qw_tot_aw
  e_val = (0.7859 + 0.03477*(t_surf-273.15))/(1.0 + 0.00412*(t_surf-273.15))
  e_sat = 10.0**e_val
  e_rat = e_sat/p_surf

  r_sat = e_rat*0.62197/(1.0 - e_rat)
  q_sat = r_sat/(1.0 + r_sat)

  !FIND MIXING RATION AT 10M HEIGHT
  r_10m = q_rel*r_sat
  q_10m = r_10m/(1.0 + r_10m)

  !GET ATM AND OCEAN HEAT FLUXES WITHOUT WIND
  q_lat_a   = rho_a*c_e*rlv*(q_sat - q_10m)
  q_sen_a   = rho_a*c_h*(t_surf - t_10m)*cp_a
  qw_tot_aw = (q_lat_a + q_sen_a)/(rho_w*cp_w)

  IF(l_root) WRITE(6,7676) qw_tot_aw !l_root = .FALSE. defined in main/les_mpi unless it is the first process

  dtdzf(1)=0.010
  dtjump  = 0.
  divgls  = 0.
  zo      = 0.0001
  zi      = -ihb !setting mixing layer
  xl      = 320.0
  yl      = 320.0
  zl      = -96.0
  izi     = NINT((zi/zl)*nnz) !finding the closest integer via rounding for the mixing layer depth

  !IF STRETCHED GRID SPECIFY LOCATION OF FIRST POINT
  zw1 = -0.5

  !DONELAN PARAMETERS (only used in skoers_ker and stokesv)
  ann      = 0.00615
  bnn      = 1.0
  f2w      = 0.13
  f_p      = f2w*grav/u_10
  npm      = 4
  sigma_p  = 2.0*pi*f_p

  time  = 0.0

  !OUTERMOST COARSE GRID INDICES ARE BOUNDS OF GRID
  izlow = 1
  izup  = nnz
  dz    = zl/nnz
  dzg   = ABS(dz)

  IF(l_root) WRITE(6,4040) zl,nnz,dzg

  !GENERATE Z GRIDS FOR PARTICULAR MESH FROM IZ = 0,..,NNZ+1
  !THIS ALLOWS INDEXING TO ARRAY ELEMENTS Z(0), ETC.
  zwstrt = 0.0

  !IF UNIFORM VERTICAL SPACING
  IF(iz_space == 0) THEN

    !BUILD Z GRID FOR W POINTS
    DO iz=0,nnz+1
      z(iz) = dz*FLOAT(iz) + zwstrt
    ENDDO
  ELSE
    CALL vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
  ENDIF

  CALL get_dz

  IF(l_root) THEN
    WRITE(6,8002) zwstrt
    WRITE(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
  ENDIF

  nnzm1 = nnz-1
  dx    = xl/nnx
  dy    = yl/nny
  fnxy  = 1./FLOAT(nxy)
  dzdz  = dzw(1)*dzw(1)
  z1    = zz(1)

  c23  = 2.0/3.0
  dsl  = (dx*1.5*dy*1.5*ABS(dzw(1)))**(1./3.)
  dslg = dsl
  cs   = 0.2

  vk     = 0.4 !von karman constant
  batagk = batag*vk
  vkin   = 1./vk
  ttmean = 0.
  zody   = ALOG(ABS(z1/zo))

  WRITE(nprt, 9901) z1,zo,zody

  zodyin = 1./zody
  wstar  = ABS(batag*zi*wtsfc(1))**(1./3.)
  !if statement below is outdated (probably from another research paper)
  IF(ismlt == 1) THEN !ismlt = 1 ; use businger formulas in MO
    !SET CONSTANTS FOR BUSINGER SIMILARITY FUNCTIONS
    vk74   = vk*0.74
    vk74in = 0.74/vk !in =inverse
    zody74 = zody*0.74
  ELSE !=0 ; use large and everyone elses formulas in MO
    !SET CONSTANTS FOR LARGE SIMILARITY FUNCTIONS
    vk74    = vk
    vk74in  = 1.0/vk !inverse
    zody74  = zody
  ENDIF

  ugal   = ugcont*0.5
  cdbtm  = vk*vk/zody/zody

  !SET SURFACE FRICTION VELOCITY HERE AND IN SR. SUFTO
  utau = SQRT(rho_a*(8.5e-4)*ws10*ws10/rho_w) !friction velocity=sqrt(tau_w/rho). tau_w=rho_a*ws10*ws10 (units check). ws10= wind speed. From Peter's turbulence textbook. 8.5e-4 is proably a drag coefficient

  utau2    = utau*utau
  IF(ibuoy == 0 .or. qstar(1) == 0.) THEN
    amonin = 1000.0
  ELSE
    amonin = -utau2*utau/(batagk*qstar(1))
  ENDIF

  hol   = ABS(zi)/amonin
  zol   = ABS(z1)/amonin
  uwsfc = -utau*utau
  vwsfc = -utau*utau

  !MAKE SURE TSFCC IS GT T00 FOR BOTH ISFC = 0,1
  tsfcc(1) = iTsurf + 273.15

  IF(l_root) THEN
    WRITE(6,80)
    WRITE(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont
  ENDIF

  IF(l_debug) THEN
    WRITE(nprt,80)
    WRITE(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont
  ENDIF

  RETURN

!FORMAT
2     FORMAT(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/, 10x,       &
            ' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,' ZO(BTM) =',e12.4, &
            ',  CDBTM=',e12.4,',  UG = ',e12.4)
80    FORMAT(///,' ***** SCRATCH RUN ***** ',//)
4040  FORMAT(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
4043  FORMAT(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
8002  FORMAT(' zwstrt = ',e12.4)
8003  FORMAT(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))

7676  FORMAT(' in init qw_tot_aw = ',e15.6)

9901  FORMAT(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,' zody = ',e15.6)

END SUBROUTINE
