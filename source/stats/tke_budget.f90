SUBROUTINE tke_budget
! GET TERMS IN RESOLVED SCALE TKE BUDGET AS IN GABLS WRITEUP AT W-POINTS AT
! ISTAGE = 1
! T_DISS AND TR_TAU ARE IN COMP1

  USE pars
  USE fields
  USE con_data
  USE con_stats

  REAL :: pxym(1:nnz), stat(1:nnz,2)

  DO iz=1,nnz
    pxym(iz)   = 0.0
    stat(iz,1) = 0.0        ! TKE TRANSPORT = WQ
    stat(iz,2) = 0.0        ! PRESSURE TRANSPORT = WP
  ENDDO

  DO iz=izs,ize

    ! GET MEAN P_STAR PRESSURE
    izm1 = iz - 1
    pxym(iz) = 0.0
    DO iy=iys,iye
      DO ix=1,nnx
        pxym(iz) = pxym(iz) + p(ix,iy,iz)-(e(ix,iy,iz)+e(ix,iy,izm1))/3.0   &
              -0.5*((u(ix,iy,iz)+stokes(iz)*dir_x)**2 + (v(ix,iy,iz)+       &
              stokes(iz)*dir_y)**2 + 0.5*(w(ix,iy,iz)*w(ix,iy,iz)+          &
              w(ix,iy,izm1)*w(ix,iy,izm1)))
      ENDDO
    ENDDO
    pxym(iz) = pxym(iz)*fnxy
  ENDDO

  CALL mpi_sum_z(pxym(1),i_root,myid,nnz,1)

  ! GET TRANSPORT TERMS AS VERTICAL ARRAYS 
  DO iz=izs,ize
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx

        ! GET ESTIMATE OF TURBULENT TRANSPORT TERM
        ufluc   = u(ix,iy,iz) - uxym(iz)
        vfluc   = v(ix,iy,iz) - vxym(iz)
        wfluc   = w(ix,iy,iz) - wxym(iz)
        wfluc_l = w(ix,iy,izm1) - wxym(izm1)
        stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*(ufluc**2 +       &
              vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))

        ! GET ESTIMATE OF PRESSURE TRANSPORT TERM
        pfluc = p(ix,iy,iz) - pxym(iz) -(e(ix,iy,iz)+e(ix,iy,izm1))/3.0 -   &
              0.5*((u(ix,iy,iz)+stokes(iz)*dir_x)**2 + (v(ix,iy,iz)+        &
              stokes(iz)*dir_y)**2 + 0.5*(w(ix,iy,iz)*w(ix,iy,iz)+          &
              w(ix,iy,izm1)*w(ix,iy,izm1)))
        stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
      ENDDO
    ENDDO
    stat(iz,1) = stat(iz,1)*fnxy
    stat(iz,2) = stat(iz,2)*fnxy
  ENDDO

  CALL mpi_sum_z(stat(1,1),i_root,myid,nnz*2,1)

  ! WE HAVE ALL TERMS ON ALL PROCESSORS FOR ALL Z, ADD THEM UP
  DO iz=1,nnz
    izp1 = iz + 1
    izm1 = iz - 1

    ! TREAT TR_TAU AT BOTTOM SPECIAL, TR_TAU = 0.0
    IF(iz == 1) tr_tau(izm1) = 0.0

    IF(iz == nnz) THEN
      t_tau(iz) = 0.0
      t_wp(iz)  = 0.0
      t_wq(iz)  = 0.0
    ELSE
      t_tau_u   = 0.5*(tr_tau(izp1) + tr_tau(iz))
      t_tau_l   = 0.5*(tr_tau(izm1) + tr_tau(iz))
      t_tau(iz) = -(t_tau_u - t_tau_l)*dzu_i(izp1)
      t_wq(iz)  = -(stat(izp1,1) - stat(iz,1))*dzu_i(izp1)
      t_wp(iz)  = -(stat(izp1,2) - stat(iz,2))*dzu_i(izp1)
    ENDIF

    dudz = (uxym(izp1) - uxym(iz))*dzu_i(izp1)
    dvdz = (vxym(izp1) - vxym(iz))*dzu_i(izp1)

    ! GATHER ALL BUDGET TERMS
    t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
    t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
    t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
    t_buoy(iz)  =  batag*wtle(iz,1)
  ENDDO

  RETURN
END SUBROUTINE
