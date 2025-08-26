SUBROUTINE recv_ytox(f, ft, nx, ixs, ixe, iys, iye, izs, ize)

    REAL :: f(nx, iys:iye, izs:ize), ft(iys:iye, ixs:ixe, izs:ize)

    DO k = izs, ize
        DO i = ixs, ixe
            DO j = iys, iye
                f(i, j, k) = ft(j, i, k)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
