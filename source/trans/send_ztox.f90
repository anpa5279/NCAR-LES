SUBROUTINE send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
!GRAB CORRECT CHUNK OF ARRAY TO BE SENT, ACCOUNTING FOR GHOST POINTS

  REAL :: g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)

  DO j=iys,iye
    DO i=ixs,ixe
      DO k=izs-1,ize+1
        gt(k,j,i) = g(k,j,i)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
