SUBROUTINE sufto(it)

  USE pars
  USE inputs
  USE fields
  USE con_data
  USE con_stats

  REAL :: buf(3+nscl)

  ! VERSION OF SIMILARITY THEORY ADAPTED FOR OCEAN FLOWS
  ! OPTION TO USE BUSINGER OR LARGE VERSION OF SIMILARITY THEORY
  iz    = 1
  izm1  = iz - 1
  izp1  = iz + 1
  z1_a  = ABS(z1)
  buf(1)  = 0.0
  buf(2)  = 0.0
  buf(3)  = 0.0
  tol     = 0.01
  DO iy=iys,iye
    DO ix=1,nnx
      buf(1) = buf(1) + u(ix,iy,iz)
      buf(2) = buf(2) + v(ix,iy,iz)
      wind(ix,iy) = SQRT((u(ix,iy,iz)+ugal)**2+v(ix,iy,iz)*v(ix,iy,iz))
      buf(3) = buf(3) + wind(ix,iy)
    ENDDO
  ENDDO

  DO iscl=1,nscl
    buf(3+iscl) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
      ENDDO
    ENDDO
  ENDDO

  ! GET X-Y SLAB SUMS
  CALL mpi_sum_xy(buf,myid,iss,ise,(3+nscl))

  u1xy  = buf(1)*fnxy + ugal
  v1xy  = buf(2)*fnxy
  windm = buf(3)*fnxy
  DO iscl=1,nscl
    t1xy(iscl) = buf(3+iscl)*fnxy
  ENDDO

  vsfc  = SQRT(u1xy*u1xy+v1xy*v1xy)
  windm = AMAX1(windm,ufree)
  vsfc  = AMAX1(vsfc,ufree)

  t10xy(1)=-qstar(1)/utau*zody*vk74in

  ! CHECK FOR TEMPERATURE BC
  IF(isfc == 0 ) THEN
    tsfcc(1)=t1xy(1)-t10xy(1)
  ENDIF

  ! INPUT SURFACE WIND STRESS (TAU = 0.0184N/M^2)
  utau = stokes(iz)*turb_la**2!sqrt(rho_a*(8.5e-4)*ws10*ws10/rho_w)

  ! SAVE OLD TAU
  utausv = utau
  utau2  = utau*utau
  IF (ibuoy==0 .OR. qstar(1) == 0.) THEN
    amonin    = 1000.
    zeta      = 0.
    thstar(1) = 0.0
    t10xy(1)  = 0.0
  ELSE
    amonin = -utau2*utau/(batagk*qstar(1))
    zeta   = z1_a/amonin
  ENDIF

  IF (t10xy(1)<0. .AND. qstar(1) .lt. 0.) THEN
    WRITE(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it

    ! TROUBLE IN SL ROUTINE
    WRITE(nprt,9000)
    CALL mpi_finalize(ierr)
    STOP
  ELSE
    ! FOR STABLE, NEUTRAL, AND UNSTABLE PBL, GET DRIFT VELOCITY
    IF(ismlt == 1) THEN
      CALL busngr(zeta,phim,phis,psim,psis)
    ELSE
      CALL fzol(zeta,phim,phis,psim,psis)
    ENDIF

    udrIFt = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
    vdrIFt = 0.0
    dnom      = (zody-psis)*vk74in
    IF (isfc==1) THEN
      thstar(1) = (t1xy(1) - tsfcc(1))/dnom
      t10xy(1)  = thstar(1)*dnom
      qstar(1)  = - utau*thstar(1)
    ELSE
      thstar(1)  = -qstar(1)/utau
      tsfcc(1)   = t1xy(1)-thstar(1)*dnom
      t10xy(1)   = thstar(1)*dnom
    ENDIF

    zol = zeta
    hol = zol*zi/z1

    ! EXAMPLES OF TWO OTHER SCALARS
    ! NOTE ROUNDOFF PROBLEM IN ANGLES ARE CLOSE TO MULTIPLES OF PI
    utau2 = utau*utau
    au13m = utau2
    au23m = 0.0
    aut3m(1)= wtsfc(1)
  ENDIF

  RETURN

! FORMAT
9000  FORMAT(' Trouble in SR. sufto')
1234  FORMAT(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,' l=',    &
            e12.3,' u*=',e12.3,' at it=',i5)

END SUBROUTINE
