SUBROUTINE setup(it)

  USE pars
  USE fields
  USE fftwk
  USE con_DATA
  USE con_stats

  it = iti
  it_counter = it - iti

  ! TURN ON NEW SGS MODEL AT A PARTICULAR STEP
  IF(it >= new_vis .AND. ivis0 == 1) THEN
    ivis = 1
  ELSE
    ivis = 0
  ENDIF

  IF(igrdr == 3) THEN
    IF(l_root) THEN
      WRITE(6,6)iti,utau,tsfcc(1) ,qstar(1)
      WRITE(6,510)
      WRITE(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont

      CALL print(6,it,1,nnz)
    ENDIF

    IF(l_debug) THEN
      WRITE(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
      WRITE(nprt,510)
      WRITE(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo,cdbtm,ugcont

      CALL print(nprt,it,izs,ize)
    ENDIF
  ENDIF

  IF(l_root) THEN
    WRITE(6,1) nnx,nny,nnz,ismlt,iti,itmax,iupwnd,ibuoy,itcut,dt,zo,        &
          tsfcc(1),method, ivis
  ENDIF

  IF(l_debug) THEN
    WRITE(nprt,1) nnx,nny,nnz,ismlt,iti,itmax,iupwnd,ibuoy,itcut,dt,zo,     &
          tsfcc(1),method, ivis
  ENDIF

  ! BC FLAGS
  ibcu = iradup
  ibcl = 0

  ! WAVENUMBERS, INTRODUCE A NORMALIZED SET OF WAVENUMBERS TO ELIMINATE
  ! COMPUTATION IN DERIVATIVES XDERIV AND YDERIV
  DO i=1,nnx
    xkn(i) = FLOAT(i-1)*pi2/xl
    IF(i>ncx)xkn(i) = -FLOAT(nnx-i+1)*pi2/xl
  ENDDO

  fn = 1.0/FLOAT(nnx)
  DO i=1,nnx
    xk(i) = xkn(i)*fn
  ENDDO

  DO i=1,nny
    ykn(i) = FLOAT(i-1)*pi2/yl
    IF(i>ncy)ykn(i) = -FLOAT(nny-i+1)*pi2/yl
  ENDDO

  fn = 1.0/FLOAT(nny)
  DO i=1,nny
    yk(i) = ykn(i)*fn
  ENDDO

  ii = -1
  DO i=1,ncx
    ii = ii + 2
    temp = xkn(i)**2
    DO j=1,nny
      temp1       = temp + ykn(j)**2
      xks(ii,j)   = temp1
      xks(ii+1,j) = temp1
    ENDDO
  ENDDO

  xnn = ABS(batag*dtdzf(1))

  ! CHOOSE CORRECT SIGN SO GRAVITY WAVES PROPAGATE OUT OF THE DOMAIN
  sgn = -1.0
  IF(ibcu==1) THEN
    DO iy=1,nny
      DO ix=1,nnxp2
        IF(xks(ix,iy) .le. 0.) THEN
          wavexy(ix,iy) = 0.0
        ELSE
          wavexy(ix,iy) = sgn*SQRT(xnn/xks(ix,iy))
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! SET LENGTH SCALE FOR SGS MODEL
  IF(iz_space == 0) THEN

    ! UNIFORM VERITCAL SPACING
    dx32 = dx*3./2.
    dy32 = dy*3./2.
    dsl  = (ABS(dx32*dy32*dzw(1)))**(1./3.)
    dslg = dsl

    IF(l_root)  WRITE(6,2000) dsl
    IF(l_debug) WRITE(nprt,2000) dsl

    ! CREATE DSL ARRAY FOR EASY INDEXING IN COMP1
    DO iz=0,nnzp1
      dsl_z(iz) = dslg
    ENDDO

  ! VARIABLE VERTICAL SPACING
  ELSE

    ! JUST ESTIMATE DSL FOR AVERAGE SPACING
    dx32 = dx*3./2.
    dy32 = dy*3./2.

    dsl_max = (ABS(dx32*dy32*dzw(0)))**(1./3.)
    DO iz=0,nnzp1
      dsl_z(iz) = (ABS(dx32*dy32*dzw(iz)))**(1./3.)
      IF(dsl_z(iz) > dsl_max) dsl_max = dsl_z(iz)
    ENDDO

    dsl  = dsl_max
    dslg = dsl
  ENDIF

  gridr = 1.0
  sml_eg = smal_e*gridr

  ! GET VISCOSITY MODEL PARAMETERS
  IF(ivis /= 1) THEN
    viscon = 0.0
    xksurf = 0.0
    nmatch = -1
    myid_newvis = 0
    DO iz=1,nnz
      dfac(iz) = 1.0
    ENDDO
  ENDIF

  ! SET STOKES VELOCITY FOR OCEANIC FLOW
  CALL stokesv
  !  CALL stokes_ideal

  ! CAN ADD A TIME SO AS TO SKIP INTO ANY PART OF THE SPECIFIED GEOSTROPIC ARRAYS
  ! ARRAYS
  ! TIME FACTOR IN SECONDS
  t_factor = 7200.0

  ! FOR PRINT OUT TO GET MORE DIGITS
  t_ref = 273.15 !+ iTsurf
  print*, ' t_ref = ', t_ref

  ! DO NOT LOOK FOR ZI BELOW ZI_MIN
  zi_min = -5.0
  iz_min = 1
  DO iz=1,nnz-1
    IF(zz(iz) < zi_min .AND. zz(iz+1) >= zi_min) iz_min = iz
  ENDDO

  IF(l_root) THEN
    WRITE(6,9000) zi_min, iz_min
  ENDIF

  RETURN

! FORMAT STATEMENTS
6 FORMAT(///,' DATA FROM RESTART FILE AT STEP =',I5,' U_* = ',e15.6,        &
      ' TS = ',e15.6,' Q_* = ',e15.6,///)
510 FORMAT(' RESTART ***** CASE WITH : ******',/)
520 FORMAT(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,', DTDZ FREE = ', &
      e12.4,', ZODY = ',e12.4,/,10x,'  ZO(BTM) = ',e12.4,', CDBTM = ',      &
      e12.4,', UG = ',e12.4)
1 FORMAT(10x,' NNX = ',i5,',  NNY = ',i5,',  NNZ = ',i5,/,10x,              &
      ' SFC SMLT = ',i1,',  ITI = ',i6,',  ITMAX = ',i6,/,10x,              &
      ' IUPWIND = ',i1,',  BUYNCY = ',i1,',  ITCUT = ',i1,/,10x,            &
      ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,10x,',  METHOD = ',  &
      i1,',  IVIS = ',i1)
2000 FORMAT(10x,' DSL = ',e15.6)
9000 FORMAT(' Search for zi above the height = ',e15.6,/,' iz_min = ',i5)

END SUBROUTINE
