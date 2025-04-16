SUBROUTINE recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)

  REAL :: g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)

  DO k=izs,ize
    DO j=iys,iye
      DO i=ixs,ixe
        g(k,j,i) = gt(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
