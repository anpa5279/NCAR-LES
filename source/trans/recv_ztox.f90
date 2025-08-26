SUBROUTINE recv_ztox(f, ft, nx, ixs, ixe, iys, iye, izs, ize)

    REAL :: f(nx, iys:iye, izs - 1:ize + 1), ft(izs - 1:ize + 1, iys:iye, ixs:ixe)

    DO i = ixs, ixe
        DO j = iys, iye
            DO k = izs - 1, ize + 1
                f(i, j, k) = ft(k, j, i)
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE
