SUBROUTINE surfvis(it)
!used if the vortex position has been updated (called in les_mpi)
!surface viscosity
  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  INCLUDE 'mpif.h'

  REAL :: xkvis(nnx,iys:iye), alwk(nnx,iys:iye)
  REAL :: send(3), buf(3)

  xksurf = 0.0
  viscon = 0.0
  vise   = 0.0

  ! ONLY ROOT PROCESS(ES) COMPUTE
  IF(iss == 0) THEN
    iz   = 1
    izm1 = iz - 1
    izp1 = iz + 1
    dz_i = dzu_i(izp1)

    CALL sufto(it)

    IF(qstar(1) == 0.) THEN !if hflux=0 
      zeta = 0.0
    ELSE
      zeta = ABS(z(1))/amonin 
    ENDIF

    IF(ismlt == 1) THEN !ismlt = 1 ; use businger formulas in MO
      CALL busngr(zeta,phim,phis,psim,psis)
    ELSE ! = 0 ; use large and everyone elses formulas in MO
      CALL fzol(zeta,phim,phis,psim,psis)
    ENDIF

    viscon = vk*ABS(z(1))/(utau*phim)
    vise   = utau*vk*ABS(z(1))/phim

    ! GET SPECIAL VALUE AT Z1 TO MATCH WITH SURFACE LAYER
    uws = 0.0
    vws = 0.0

    DO iy=iys,iye
      DO ix=1,nnx
        uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) + u(ix,iy,izp1) - u_mn(izp1)) &
              *(w(ix,iy,iz)-w_mn(iz))
        vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) + v(ix,iy,izp1) - v_mn(izp1)) &
              *(w(ix,iy,iz)-w_mn(iz))
      ENDDO
    ENDDO

    uws = uws*fnxy
    vws = vws*fnxy

    ! GET AVERAGE FLUCTUATING EDDY VISCOSITY
    DO iy=iys,iye
      DO ix=1,nnx
        e(ix,iy,iz)=AMAX1(e(ix,iy,iz),sml_eg)
      ENDDO
    ENDDO

    dslk = AMIN1(dsl,vk*ABS(z(iz))/csmag)
    almin = almin_c*dsl_z(iz)
    DO iy=iys,iye
      DO ix=1,nnx
        alwk(ix,iy)=dslk

        ! NO STABILITY CORRECTED LENGTH SCALES WHEN NEW EDDY VISCOSITY IS ON
        xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
      ENDDO
    ENDDO

    ! GET AVERAGE VISCOSITY
    xkavg = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        xkavg = xkavg + xkvis(ix,iy)
      ENDDO
    ENDDO

    xkavg = xkavg*fnxy

    buf(1) = uws
    buf(2) = vws
    buf(3) = xkavg

    CALL mpi_sum_xy(buf,myid,iss,ise,3)

    uws   = buf(1)
    vws   = buf(2)
    xkavg = buf(3)

    xkz1 = vise - SQRT(uws**2 + vws**2)*viscon
    xksurf =  xkz1 - xkavg
    xksurf = AMAX1(xksurf,0.0) !ensuring it is not negative. 
    xksurf = AMIN1(xksurf,vise) !can be postiive or negative
  ENDIF

  ! BROADCAST VALUES TO OTHER PROCESSES
  send(1) = xksurf
  send(2) = viscon
  send(3) = vise

  CALL mpi_bcast(send,3,mpi_real8,i_root,mpi_comm_world,ierr)

  xksurf = send(1)
  viscon = send(2)
  vise   = send(3)

  RETURN

! FORMAT
6000  FORMAT(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,         &
            ' vise = ',e12.4,' xksurf = ',e12.4)

END SUBROUTINE
