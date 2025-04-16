SUBROUTINE tridv(b,d,a,r,n,j1,j2)!tridv(bb,dd,aa,rh,nnz,jxs,jxe)
!TRIDIAGONAL MATRIX SOLVER WITH MULTIPLE VECTORS (NOTE: J AND I LOOPS ARE
!REVERSED FROM CRAY VERSION)
!INPUT:  N     - SIZE OF A,B,D, AND R (nnz)
!        B     - BELOW DIAGONAL ELEMENTS (B(1) NOT USED) (bb)
!        D     - DIAGONAL ELEMENTS (dd)
!        A     - ABOVE DIAGONAL ELEMENTS (A(N) NOT USED) (aa)
!        R     - RHS (rh)
!        J1:J2 - RANGE OF INPUT VECTORS (jxs, jxe)
!OUTPUT: R     - SOLUTION VECTOR (rh)
!only called in solve/solve_trid tridv(bb,dd,aa,rh,nnz,jxs,jxe)

  REAL :: b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)

  IF(n > 1) THEN
    DO j=j1,j2
      d(1,j) = 1.0/d(1,j)
    ENDDO

    DO j=j1,j2
      DO i=2,n
        fac = b(i,j)*d(i-1,j)
        d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
        r(i,j) = r(i,j) - fac*r(i-1,j)
      ENDDO
    ENDDO

    DO j=j1,j2
      r(n,j) = r(n,j)*d(n,j)
    ENDDO

    DO j=j1,j2
      DO i=n-1,1,-1
        r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      ENDDO
    ENDDO
  ELSE
    DO j=j1,j2
      r(1,j) = r(1,j)/d(1,j)
    ENDDO
  ENDIF

  RETURN
END SUBROUTINE
