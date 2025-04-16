MODULE fields

  IMPLICIT NONE

  REAL, ALLOCATABLE ::                                                      & !velocity, energy, and r for scalar
          u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:), r1(:,:,:),    &
          r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:)
  REAL, ALLOCATABLE ::                                                      &
          ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:), wx(:,:,:), wy(:,:,:), & !spatial derivatives of velocity, pressure (what is ptop?), viscosity
          p(:,:,:), ptop(:,:,:), vis_m(:,:,:), vis_s(:,:,:), vis_sv(:,:,:)
  REAL, ALLOCATABLE ::                                                      &
          ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:), ebc(:,:,:),     & !velocity, scalar, energy, and pressure boundary conditions (updates the first set of 3D matrices above)
          pbc(:,:,:), pbc2(:,:,:)

END MODULE
