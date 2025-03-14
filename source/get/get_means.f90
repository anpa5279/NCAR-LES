SUBROUTINE get_means(istage)
! GET MEANS FOR ALL VARIABLES FOR USE IN ISO, SURFVIS, COMP1, COMPMN

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  DO iz=0,nnz+1
    u_mn(iz)   = 0.0
    v_mn(iz)   = 0.0
    w_mn(iz)   = 0.0
    engz(iz)   = 0.0
    engsbz(iz) = 0.0
    divz(iz)   = 0.0
  ENDDO

  DO iscl=1,nscl
    DO iz=0,nnz+1
      t_mn(iz,iscl) = 0.0
    ENDDO
  ENDDO

  iz_ee = ize
  IF(ize == nnz) iz_ee = nnzp1

  DO iz=izs,iz_ee
    DO iy=iys,iye
      DO ix=1,nnx
        u_mn(iz) = u_mn(iz) + u(ix,iy,iz)
        v_mn(iz) = v_mn(iz) + v(ix,iy,iz)
        w_mn(iz) = w_mn(iz) + w(ix,iy,iz)
      ENDDO
    ENDDO

    u_mn(iz) = u_mn(iz)*fnxy
    v_mn(iz) = v_mn(iz)*fnxy
    w_mn(iz) = w_mn(iz)*fnxy
    DO iscl=1,nscl
      t_mn(iz,iscl) = 0.0
      DO iy=iys,iye
        DO ix=1,nnx
          t_mn(iz,iscl) = t_mn(iz,iscl) + t(ix,iy,iscl,iz)
        ENDDO
      ENDDO

      t_mn(iz,iscl) = t_mn(iz,iscl)*fnxy
    ENDDO
  ENDDO

  CALL mpi_sum_z(u_mn(1),i_root,myid,nnzp1,1)
  CALL mpi_sum_z(v_mn(1),i_root,myid,nnzp1,1)
  CALL mpi_sum_z(w_mn(1),i_root,myid,nnzp1,1)

  DO iscl=1,nscl
    CALL mpi_sum_z(t_mn(1,iscl),i_root,myid,nnzp1,1)
  ENDDO

  ! SET E TO MINIMUM VALUE
  DO iz=izs-1,ize+1
    DO iy=iys,iye
      DO ix=1,nnx
        e(ix,iy,iz)=AMAX1(e(ix,iy,iz ),sml_eg)
      ENDDO
    ENDDO
  ENDDO

  ! GET TERMS WHICH CONTRIBUTE TO MEAN PRESSURE CAREFUL WITH SUM
  DO iz=izs,ize
    izm1 = iz - 1
    DO iy=iys,iye
      DO ix=1,nnx
        engz(iz)  = engz(iz) + .5*((u(ix,iy,iz) + stokes(iz)*dir_x)**2 +  &
              (v(ix,iy,iz) + stokes(iz)*dir_y)**2 + .5*(w(ix,iy,iz)*      &
              w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
        engsbz(iz) = engsbz(iz) + 0.5*(e(ix,iy,iz)+e(ix,iy,izm1))
      ENDDO
    ENDDO
    engz(iz)   = engz(iz)*fnxy
    engsbz(iz) = engsbz(iz)*fnxy
  ENDDO

  CALL mpi_sum_z(engz(1),i_root,myid,nnzp1,1)
  CALL mpi_sum_z(engsbz(1),i_root,myid,nnzp1,1)

  ! SAVE MEANS AND DIVERGENCE FOR PRINTOUT AND COMPMN
  ! ALL CPUS HAVE MEANS OVER ALL Z
  IF(istage .eq. 1) THEN
    DO iz=izs,ize
      izm1 = iz - 1
      DO iy=iys,iye
        DO ix=1,nnx
          divz(iz) = divz(iz) + (ux(ix,iy,iz)+vy(ix,iy,iz)+(w(ix,iy,iz)-    &
                w(ix,iy,izm1))*dzw_i(iz))**2
        ENDDO
      ENDDO
      divz(iz) = divz(iz)*fnxy
    ENDDO

    CALL mpi_sum_z(divz(1),i_root,myid,nnz,1)

    DO iz=0,nnz+1
      uxym(iz) = u_mn(iz)
      vxym(iz) = v_mn(iz)
      wxym(iz) = w_mn(iz)
    ENDDO

    DO iscl=1,nscl
      DO iz=0,nnz+1
        txym(iz,iscl) = t_mn(iz,iscl)
      ENDDO
    ENDDO
  ENDIF

  RETURN
END SUBROUTINE
