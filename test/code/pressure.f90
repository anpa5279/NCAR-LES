SUBROUTINE pressure
!SOLVE FOR PRESSURE USING A MATRIX TRANSPOSE ACROSS MPI TASKS AND TRIDIAG
!SOLVER.
!THE TRANSPOSED ARRAY IS DIMENSIONED (0:NNZ+1).
!VALUES (0 & NNZ+1) ARE NOT NEEDED BUT ARE USEFUL IN THE MATRIX TRANSPOSE WHEN
!WE RETURN (SEE SEND_ZTOX).
!ON EXIT P IS DEFINED AT ALL [IZS-1:IZE+1]
!only called in solve/comp_p

  USE pars
  USE fields
  USE fftwk
  USE con_data
  USE con_stats

  REAL :: pfft(nny,jxs:jxe,izs-1:ize+1)
  REAL :: pt(0:nnz+1,jxs:jxe,jys:jye)
  REAL :: ptopfft(nny,jxs:jxe,1:2)
  REAL :: psum(1:nnz)

  INCLUDE 'mpif.h'

  INTEGER :: istatus(mpi_status_size)

  !FOURIER ANALYZE THE RHS AT ALL IZ = IZS,IZE.
  !RESULTS ARE IN PFFT
  CALL fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,nnx,nny,jxs, &
        jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs,-2)

  !FOURIER ANALYZE THE RADIATION BC ARRAYS
  IF(ibcu == 1) THEN
    CALL fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),trigx(1,1),trigc,nnx,nny, &
          jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,1,2,myid,ncpu_s,numprocs,-2)
  ENDIF

  !TRANSPOSE FIRST AND LAST INDEX OF ARRAY
  !THE ORDER OF PFFT IS (Y,X,Z)
  CALL xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,jxs,jxe,izs,ize,iz_s,   &
        iz_e,myid,ncpu_s,numprocs)
  CALL solve_trid(pt, ptopfft) !only time this is called

  !TRANSPOSE BACK
  CALL ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,jxs,jxe,izs,ize,iz_s,   &
        iz_e,myid,ncpu_s,numprocs)

  iz_ee = ize+1
  IF(ise == numprocs-1) THEN
    iz_ee = ize
  ENDIF

  !INVERSE FFT AT ALL IZ = IZS, IZ_EE TO GET P
  !SEE Z INDICES
  CALL fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,nnx,nny,jxs, &
        jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,izs,iz_ee,myid,ncpu_s,numprocs,2)

  !PARTIAL SUMS FOR PRESSURE
  DO iz=1,nnz
    psum(iz) = 0.0
  ENDDO

  DO iz=izs,ize
    DO iy=iys,iye
      DO ix=1,nnx
        psum(iz) = psum(iz) + p(ix,iy,iz)
      ENDDO
    ENDDO
    psum(iz) = psum(iz)*fnxy
  ENDDO

  CALL mpi_sum_z(psum,i_root,myid,nnz,1)

  DO iz=izs,iz_ee
    DO iy=iys,iye
      DO ix=1,nnx
        p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE
