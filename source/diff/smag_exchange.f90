subroutine smag_exchange
! exchange ghost points with mpi,
! nb and nt are the destination and
! source nodes. Allows for 1z per cpu

    use pars
    use fields

    include 'mpif.h'

    real fs(nnx, iys:iye), fr(nnx, iys:iye)
    integer istatus(MPI_status_size)

    nb = myid - ncpu_s
    nt = myid + ncpu_s

    ! account for endpoints
    if (iss == 0) then
        nb = MPI_proc_null
    end if
    if (ise == numprocs - 1) then
        nt = MPI_proc_null
    end if
    nsend = nnx * (iye - iys + 1)
    nrecv = nsend

    ! send top of myid, receive bottom from myid - ncpu_s
    do iy = iys, iye
        do ix = 1, nnx
            fs(ix, iy) = vis_sv(ix, iy, ize)
        end do
    end do

    call mpi_sendrecv(fs(1, iys), nsend, MPI_real8, nt, 0, &
                      fr(1, iys), nrecv, MPI_real8, nb, 0, MPI_comm_world, istatus, ierr)

    if (iss /= 0) then
        do iy = iys, iye
            do ix = 1, nnx
                vis_sv(ix, iy, izs - 1) = fr(ix, iy)
            end do
        end do
    end if

! send bottom of myid, receive bottom from myid + ncpu_s
    do iy = iys, iye
        do ix = 1, nnx
            fs(ix, iy) = vis_sv(ix, iy, izs)
        end do
    end do

    call mpi_sendrecv(fs(1, iys), nsend, MPI_real8, nb, 1, &
                      fr(1, iys), nrecv, MPI_real8, nt, 1, MPI_comm_world, istatus, ierr)

    if (ise /= numprocs - 1) then
        do iy = iys, iye
            do ix = 1, nnx
                vis_sv(ix, iy, ize + 1) = fr(ix, iy)
            end do
        end do
    end if

    return
end
