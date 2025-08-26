SUBROUTINE tridv(b, d, a, r, n, j1, j2)
! TRIDIAGONAL MATRIX SOLVER WITH MULTIPLE VECTORS (NOTE: J AND I LOOPS ARE
! REVERSED FROM CRAY VERSION)
! INPUT:  N     - SIZE OF A,B,D, AND R
!         B     - BELOW DIAGONAL ELEMENTS (B(1) NOT USED)
!         D     - DIAGONAL ELEMENTS
!         A     - ABOVE DIAGONAL ELEMENTS (A(N) NOT USED)
!         R     - RHS
!         J1:J2 - RANGE OF INPUT VECTORS
! OUTPUT: R     - SOLUTION VECTOR

    REAL :: b(n, j1:j2), d(n, j1:j2), a(n, j1:j2), r(n, j1:j2)

    IF (n > 1) THEN
        DO j = j1, j2
            d(1, j) = 1.0 / d(1, j)
        END DO

        DO j = j1, j2
            DO i = 2, n
                fac = b(i, j) * d(i - 1, j)
                d(i, j) = 1.0 / (d(i, j) - fac * a(i - 1, j))
                r(i, j) = r(i, j) - fac * r(i - 1, j)
            END DO
        END DO

        DO j = j1, j2
            r(n, j) = r(n, j) * d(n, j)
        END DO

        DO j = j1, j2
            DO i = n - 1, 1, -1
                r(i, j) = d(i, j) * (r(i, j) - a(i, j) * r(i + 1, j))
            END DO
        END DO
    ELSE
        DO j = j1, j2
            r(1, j) = r(1, j) / d(1, j)
        END DO
    END IF

    RETURN
END SUBROUTINE
