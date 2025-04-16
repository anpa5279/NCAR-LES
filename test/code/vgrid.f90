SUBROUTINE vgrid(z1,zi,zl,nnz,z,l_root,ldebug)

  REAL :: z(0:nnz+1)
  LOGICAL :: l_root, l_debug

  !BUILD GRID UP TO ZI FIRST
  z_frst = z1
  z_cntr = zl
  n_pbl  = nnz
  z_fac1 = z_cntr/z_frst
  z_fac2 = 1.0/FLOAT(nnz)
  z_fac  = 1.1
  knt = 0
  tol = 1.0e-10

  DO WHILE(test > tol)
    knt = knt + 1
    z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
    test   = ABS(1.0 - z_facn/z_fac)
    IF(knt > 500) THEN
      IF(l_root) WRITE(6,9000) z_fac, z_facn, knt
      STOP
    ENDIF
    z_fac = z_facn
  ENDDO

  IF(l_root) WRITE(6,9100) z_fac, z_cntr, z1, knt

  z(1) = z_frst
  DO iz=2,n_pbl
    z(iz) = z_frst*(z_fac**(FLOAT(iz)) - 1.0)/(z_fac - 1.0)
  ENDDO

  z(nnz) = zl
  z(0)   = 0.0
  z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

  IF(l_root) WRITE(6,5300) n_pbl

  RETURN

!FORMAT
9000  FORMAT(' Cannot find stretching factor',/,' z_fac = ',e15.6,          &
            ' z_facn = ',e15.6,' knt = ',i3)
9100  FORMAT(' Stretching factor = ',e15.6,/,' Match point       = ',       &
            e15.6,/,' First z           = ',e15.6,/,' Number of iters   = ',i4)
5300  FORMAT(' n_pbl = ',i4)
5600  FORMAT(' 5600 in vgrid ',/,' iz ',5x,' zw ',/,(i3,e15.6))

END SUBROUTINE
