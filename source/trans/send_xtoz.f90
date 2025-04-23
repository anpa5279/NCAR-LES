SUBROUTINE send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
! GRAB CORRECT CHUNK OF ARRAY TO BE SENT AND SKIP GHOST POINTS

  REAL :: f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)

  DO k=izs,ize
    DO j=iys,iye
      DO i=ixs,ixe
        ft(i,j,k) = f(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
