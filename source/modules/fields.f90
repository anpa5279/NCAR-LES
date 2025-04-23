MODULE fields

  IMPLICIT NONE

  REAL, ALLOCATABLE ::                                                      &
          u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:), r1(:,:,:),    &
          r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:)
  REAL, ALLOCATABLE ::                                                      &
          ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:), wx(:,:,:), wy(:,:,:), &
          p(:,:,:), ptop(:,:,:), vis_m(:,:,:), vis_s(:,:,:), vis_sv(:,:,:)
  REAL, ALLOCATABLE ::                                                      &
          ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:), ebc(:,:,:),     &
          pbc(:,:,:), pbc2(:,:,:)

END MODULE
