MODULE con_stats

  USE pars
  IMPLICIT NONE

! ============================================================================ !

  TYPE con_s
    SEQUENCE
    REAL ::                                                                 &
            wwsb(maxnz), engz(0:maxnz1), engsbz(0:maxnz1), englez(maxnz),   &
            uxym(0:maxnz1), vxym(0:maxnz1), wxym(0:maxnz1),                 &
            txym(0:maxnz1,nscl),divz(0:maxnz1), utle(maxnz,nscl),           &
            utsb(maxnz,nscl), vtle(maxnz,nscl), vtsb(maxnz,nscl),           &
            wtle(maxnz,nscl), wtsb(maxnz,nscl), wt_tot(maxnz,nscl),         &
            z(0:maxnz1), zz(0:maxnz1), shrz(maxnz), buyz(maxnz),            &
            triz(maxnz), uwsb(maxnz), vwsb(maxnz), uwle(maxnz),             &
            vwle(maxnz), uw_tot(maxnz), vw_tot(maxnz), wcube(maxnz),        &
            wfour(maxnz), tcube(maxnz,nscl), ups(maxnz), vps(maxnz),        &
            wps(maxnz), tps(maxnz,nscl), t_rprod(maxnz), t_wq(maxnz),       &
            t_wp(maxnz), t_tau(maxnz), t_tran(maxnz), t_buoy(maxnz),        &
            t_diss(maxnz), t_sprod(maxnz)
    REAL ::                                                                 &
            xkn(maxnx), ykn(maxny), xk(maxnx),yk(maxny), xks(maxnx2,maxny), &
            wavexy(maxnx2,maxny)
    REAL ::                                                                 &
            ug(maxnz), vg(maxnz), wls(maxnz),uls(maxnx)
    REAL ::                                                                 &
            udrift,vdrift,cpou10,turb_la, stokesw,stokesa, stokess,         &
            stokes(maxnz1),f2w,ann,gravity,sigma_p,bnn,z_pt
    REAL ::                                                                 &
            dtg, dslg, dzg
  END TYPE


  TYPE(con_s), TARGET :: c_s
    REAL, POINTER ::                                                        &
            wwsb(:), engz(:), engsbz(:), englez(:), uxym(:), vxym(:),       &
            wxym(:), txym(:,:), divz(:), utle(:,:), utsb(:,:), vtle(:,:),   &
            vtsb(:,:), wtle(:,:), wtsb(:,:), wt_tot(:,:), z(:), zz(:),      &
            shrz(:), buyz(:), triz(:), uwsb(:), vwsb(:), uwle(:), vwle(:),  &
            uw_tot(:), vw_tot(:), wcube(:), wfour(:), tcube(:,:), ups(:),   &
            vps(:), wps(:), tps(:,:), t_rprod(:), t_wq(:), t_wp(:),         &
            t_tau(:), t_tran(:), t_buoy(:), t_diss(:), t_sprod(:)
    REAL, POINTER ::                                                        &
            xkn(:), ykn(:), xk(:), yk(:), xks(:,:), wavexy(:,:)
    REAL, POINTER ::                                                        &
            ug(:), vg(:), wls(:), uls(:)
    REAL, POINTER ::                                                        &
            udrift, vdrift, cpou10, turb_la, stokesw, stokesa, stokess,     &
            stokes(:),f2w,ann,gravity,sigma_p,bnn,z_pt
    REAL, POINTER ::                                                        &
            dtg, dslg, dzg

! ---------------------------------------------------------------------------- !
  CONTAINS

  SUBROUTINE fill_cs
  IMPLICIT NONE
! POINTER ASSOCIATION FOR STAT ARRAYS
! GET SIZE OF STAT ARRAYS ISIZE FOR HISTORY FILES

    isize = 0

    wwsb    => c_s%wwsb     ; isize = isize + size(wwsb)
    engz    => c_s%engz     ; isize = isize + size(engz)
    engsbz  => c_s%engsbz   ; isize = isize + size(engsbz)
    englez  => c_s%englez   ; isize = isize + size(englez)
    uxym    => c_s%uxym     ; isize = isize + size(uxym)
    vxym    => c_s%vxym     ; isize = isize + size(vxym)
    wxym    => c_s%wxym     ; isize = isize + size(wxym)
    txym    => c_s%txym     ; isize = isize + size(txym)
    divz    => c_s%divz     ; isize = isize + size(divz)
    utle    => c_s%utle     ; isize = isize + size(utle)
    utsb    => c_s%utsb     ; isize = isize + size(utsb)
    vtle    => c_s%vtle     ; isize = isize + size(vtle)
    vtsb    => c_s%vtsb     ; isize = isize + size(vtsb)
    wtle    => c_s%wtle     ; isize = isize + size(wtle)
    wtsb    => c_s%wtsb     ; isize = isize + size(wtsb)
    wt_tot  => c_s%wt_tot   ; isize = isize + size(wt_tot)
    z       => c_s%z        ; isize = isize + size(z)
    zz      => c_s%zz       ; isize = isize + size(zz)
    shrz    => c_s%shrz     ; isize = isize + size(shrz)
    buyz    => c_s%buyz     ; isize = isize + size(buyz)
    triz    => c_s%triz     ; isize = isize + size(triz)
    uwsb    => c_s%uwsb     ; isize = isize + size(uwsb)
    vwsb    => c_s%vwsb     ; isize = isize + size(vwsb)
    uwle    => c_s%uwle     ; isize = isize + size(uwle)
    vwle    => c_s%vwle     ; isize = isize + size(vwle)
    uw_tot  => c_s%uw_tot   ; isize = isize + size(uw_tot)
    vw_tot  => c_s%vw_tot   ; isize = isize + size(vw_tot)
    wcube   => c_s%wcube    ; isize = isize + size(wcube)
    wfour   => c_s%wfour    ; isize = isize + size(wfour)
    tcube   => c_s%tcube    ; isize = isize + size(tcube)
    ups     => c_s%ups      ; isize = isize + size(ups)
    vps     => c_s%vps      ; isize = isize + size(vps)
    wps     => c_s%wps      ; isize = isize + size(wps)
    tps     => c_s%tps      ; isize = isize + size(tps)
    t_rprod => c_s%t_rprod  ; isize = isize + size(t_rprod)
    t_wq    => c_s%t_wq     ; isize = isize + size(t_wq)
    t_wp    => c_s%t_wp     ; isize = isize + size(t_wp)
    t_tau   => c_s%t_tau    ; isize = isize + size(t_tau)
    t_tran  => c_s%t_tran   ; isize = isize + size(t_tran)
    t_buoy  => c_s%t_buoy   ; isize = isize + size(t_buoy)
    t_diss  => c_s%t_diss   ; isize = isize + size(t_diss)
    t_sprod => c_s%t_sprod  ; isize = isize + size(t_sprod)
    xkn     => c_s%xkn
    ykn     => c_s%ykn
    xk      => c_s%xk
    yk      => c_s%yk
    xks     => c_s%xks
    wavexy  => c_s%wavexy
    ug      => c_s%ug
    vg      => c_s%vg
    wls     => c_s%wls
    uls     => c_s%uls
    udrift  => c_s%udrift
    vdrift  => c_s%vdrift
    vdrift  => c_s%vdrift
    cpou10  => c_s%cpou10
    turb_la => c_s%turb_la
    stokesw => c_s%stokesw
    stokesa => c_s%stokesa
    stokess => c_s%stokess
    stokes  => c_s%stokes
    f2w     => c_s%f2w
    ann     => c_s%ann
    sigma_p => c_s%sigma_p
    bnn     => c_s%bnn
    z_pt    => c_s%z_pt
    dtg     => c_s%dtg
    dslg    => c_s%dslg
    dzg     => c_s%dzg
  RETURN
  END SUBROUTINE
END MODULE
