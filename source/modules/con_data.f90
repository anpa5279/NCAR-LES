MODULE con_data

  USE pars, ONLY: nscl
  IMPLICIT NONE

! ============================================================================ !

  TYPE con_d
    SEQUENCE
    REAL ::                                                                 &
            zo, vk, vkin, vk74, vk74in, grav, gcp, fcor, fcor_h, zi, two_pi,   &
            batagk, t00, batag, vgcont, ugcont, cdbtm, dtdzf(nscl), dtjump, &
            ugal, divgls, z1, utausv, xl, yl, zl, dx, dy, dz, dt, fnxy,     &
            dzdz, dsl, c23, dtgama, dtzeta, xkmax, time, zody, zody74,      &
            tsfcc(nscl), utau, qstar(nscl), wtsfc(nscl), uwsfc, vwsfc,      &
            amonin, zol, hol, smal_e, sml_eg, rho_a, rho_w, u_10, v_10,     &
            cd_10, tau_x, tau_y
  END TYPE

  TYPE(con_d), TARGET :: c_c
    REAL, POINTER ::                                                        &
            zo, vk, vkin, vk74, vk74in, grav, gcp, fcor, fcor_h, zi, two_pi,   &
            batagk, t00, batag, vgcont, ugcont, cdbtm, dtdzf(:), dtjump,    &
            ugal, divgls, z1, utausv, xl, yl, zl, dx, dy, dz, dt, fnxy,     &
            dzdz, dsl, c23, dtgama, dtzeta, xkmax, time, zody, zody74,      &
            tsfcc(:), utau, qstar(:), wtsfc(:), uwsfc, vwsfc, amonin, zol,  &
            hol, smal_e, sml_eg, rho_a, rho_w, u_10, v_10, cd_10, tau_x,    &
            tau_y

! ---------------------------------------------------------------------------- !
  CONTAINS

  SUBROUTINE fill_cc
  IMPLICIT NONE
! POINTER ASSOCIATIONS FOR CONSTANT VARIABLES

    zo     => c_c%zo
    vk     => c_c%vk
    vkin   => c_c%vkin
    vk74   => c_c%vk74
    vk74in => c_c%vk74in
    grav   => c_c%grav
    gcp    => c_c%gcp
    fcor   => c_c%fcor
    fcor_h => c_c%fcor_h
    zi     => c_c%zi
    two_pi    => c_c%two_pi
    batagk => c_c%batagk
    t00    => c_c%t00
    batag  => c_c%batag
    vgcont => c_c%vgcont
    ugcont => c_c%ugcont
    cdbtm  => c_c%cdbtm
    dtdzf  => c_c%dtdzf
    dtjump => c_c%dtjump
    ugal   => c_c%ugal
    divgls => c_c%divgls
    z1     => c_c%z1
    utausv => c_c%utausv
    xl     => c_c%xl
    yl     => c_c%yl
    zl     => c_c%zl
    dx     => c_c%dx
    dy     => c_c%dy
    dz     => c_c%dz
    dt     => c_c%dt
    fnxy   => c_c%fnxy
    dzdz   => c_c%dzdz
    dsl    => c_c%dsl
    c23    => c_c%c23
    dtgama => c_c%dtgama
    dtzeta => c_c%dtzeta
    xkmax  => c_c%xkmax
    time   => c_c%time
    zody   => c_c%zody
    zody74 => c_c%zody74
    tsfcc  => c_c%tsfcc
    utau   => c_c%utau
    qstar  => c_c%qstar
    wtsfc  => c_c%wtsfc
    uwsfc  => c_c%uwsfc
    vwsfc  => c_c%vwsfc
    amonin => c_c%amonin
    zol    => c_c%zol
    hol    => c_c%hol
    smal_e => c_c%smal_e
    sml_eg => c_c%sml_eg
    rho_a  => c_c%rho_a
    rho_w  => c_c%rho_w

    u_10   => c_c%u_10
    v_10   => c_c%v_10
    cd_10  => c_c%cd_10
    tau_x  => c_c%tau_x
    tau_y  => c_c%tau_y

  RETURN
  END SUBROUTINE
END MODULE 
