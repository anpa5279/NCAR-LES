SUBROUTINE xy_stats
! GET STATISTICS
! only called in les_mpi
  USE pars
  USE fields
  USE con_data
  USE con_stats

  ! INDICES FOR INDEXING ARRAY STAT(.,.)
  INTEGER, PARAMETER ::                   &
      js = 10,                            & ! NUMBER OF NON-SCALAR STATS (never leaves this subrou)
      ns = 5                                ! NUMBER OF SCALAR STATS
  INTEGER, PARAMETER :: nstat = js + ns*nscl
  REAL :: stat(1:nnz,nstat)

  ! USE TRICK WITH MPI REDUCE OVER ALL Z TO GET AVERAGES BY SETTING STAT ARRAY
  ! EQUAL TO 0 FOR ALL Z ON EACH PROCESS
  DO i=1,nstat
    DO iz=1,nnz
      stat(iz,i) = 0.0
    ENDDO
  ENDDO

  ! INDICES FOR SCALARS
  m1 = js
  m2 = js + nscl
  m3 = js + 2*nscl
  m4 = js + 3*nscl
  m5 = js + 4*nscl

  sgn = 1.0
  IF(iupwnd == 1) sgn = -1.0

  DO iz=izs,ize
    izp2 = iz + 2
    izp1 = iz + 1
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx
        stat(iz,1) = stat(iz,1) + (u(ix,iy,iz) - uxym(iz))**2
        stat(iz,2) = stat(iz,2) + (v(ix,iy,iz) - vxym(iz))**2
        stat(iz,3) = stat(iz,3) + (w(ix,iy,iz) - wxym(iz))**2
        stat(iz,4) = stat(iz,4) + (w(ix,iy,iz) - wxym(iz))**3
        stat(iz,5) = stat(iz,5) + (w(ix,iy,iz) - wxym(iz))**4
        stat(iz,6) = stat(iz,6) + ((w(ix,iy,iz)-wxym(iz))**2 + (0.5*        &
              (u(ix,iy,iz)-uxym(iz) + u(ix,iy,izp1)-uxym(izp1)))**2 +       &
              (0.5*(v(ix,iy,iz)-vxym(iz) + v(ix,iy,izp1)-vxym(izp1)))**2)*0.5
        stat(iz,7) = stat(iz,7) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
        stat(iz,8) = stat(iz,8) + e(ix,iy,iz)
        stat(iz,9) = stat(iz,9) + (w(ix,iy,iz)-wxym(iz))*0.5*((u(ix,iy,iz)- &
              uxym(iz))+(u(ix,iy,izp1)-uxym(izp1)))
        stat(iz,10) = stat(iz,10) + (w(ix,iy,iz)-wxym(iz))*0.5*             &
              ((v(ix,iy,iz)-vxym(iz))+(v(ix,iy,izp1)-vxym(izp1)))
      ENDDO
    ENDDO

    ! GET SCALAR RESOLVED FLUXES AND VARIANCES
    DO l=1,nscl
      IF(iupwnd /= 1 .or. iz == nnz) THEN
        DO iy=iys,iye
          DO ix=1,nnx
            stat(iz,m1+l)=stat(iz,m1+l) + (w(ix,iy,iz)-wxym(iz))*0.5*       &
                  (t(ix,iy,l,iz)-txym(iz,l) + t(ix,iy,l,izp1)-txym(izp1,l))
          ENDDO
        ENDDO
      ELSE

        ! MONOTONE FLUXES
        DO iy=iys,iye
          DO ix=1,nnx
            stat(iz,m1+l) = stat(iz,m1+l) + AMAX1(sgn*w(ix,iy,iz),0.)*      &
                  (t(ix,iy,l,iz) + rlim(t(ix,iy,l,izp1),t(ix,iy,l,iz),      &
                  t(ix,iy,l,izm1))) + AMIN1(sgn*w(ix,iy,iz),0.)*            &
                  (t(ix,iy,l,izp1) + rlim(t(ix,iy,l,iz),t(ix,iy,l,izp1),    &
                  t(ix,iy,l,izp2)))
          ENDDO
        ENDDO
      ENDIF

      stat(iz,m1+l)= sgn*stat(iz,m1+l)

      ! GET HORIZONTAL SCALAR RESOLVED FLUXES
      DO iy=iys,iye
        DO ix=1,nnx
          stat(iz,m2+l) = stat(iz,m2+l)+(u(ix,iy,iz)-uxym(iz))*            &
                (t(ix,iy,l,iz)-txym(iz,l))
          stat(iz,m3+l) = stat(iz,m3+l)+(v(ix,iy,iz)-vxym(iz))*            &
                (t(ix,iy,l,iz)-txym(iz,l))
        ENDDO
      ENDDO

      ! SCALAR VARIANCES & HIGHER MOMENTS
      DO iy=iys,iye
        DO ix=1,nnx
          stat(iz,m4+l) = stat(iz,m4+l) + (t(ix,iy,l,iz) - txym(iz,l))**2
          stat(iz,m5+l) = stat(iz,m5+l) + (t(ix,iy,l,iz) - txym(iz,l))**3
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! ADD PARTIAL SUMS AND SEND TO ALL
  CALL mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)

  ! FILL ARRAYS FOR PRINTOUT AND CONSTANT FILE
  DO iz=1,nnz
    ups(iz)    = stat(iz,1)*fnxy
    vps(iz)    = stat(iz,2)*fnxy
    wps(iz)    = stat(iz,3)*fnxy
    wcube(iz)  = stat(iz,4)*fnxy
    wfour(iz)  = stat(iz,5)*fnxy
    englez(iz) = stat(iz,6)*fnxy
    engsbz(iz) = stat(iz,7)*fnxy
    eavg(iz)   = stat(iz,8)*fnxy
    uwle(iz)   = stat(iz,9)*fnxy
    vwle(iz)   = stat(iz,10)*fnxy
    uw_tot(iz) = uwle(iz) + uwsb(iz)
    vw_tot(iz) = vwle(iz) + vwsb(iz)

    ! GET SCALAR RESOLVED FLUXES AND VARIANCES
    DO l=1,nscl
      wtle(iz,l)   = stat(iz,m1+l)*fnxy
      utle(iz,l)   = stat(iz,m2+l)*fnxy
      vtle(iz,l)   = stat(iz,m3+l)*fnxy
      tps(iz,l)    = stat(iz,m4+l)*fnxy
      tcube(iz,l)  = stat(iz,m5+l)*fnxy
      wt_tot(iz,l) = wtle(iz,l) + wtsb(iz,l)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
