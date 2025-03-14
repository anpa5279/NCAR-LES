SUBROUTINE send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
!GRAB CORRECT CHUNK OF ARRAY TO BE SENT

  REAL :: f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)

  DO k=izs,ize
    DO j=iys,iye
      DO i=ixs,ixe
        ft(i,j,k) = f(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
