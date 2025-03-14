module reaction
  use fields, only: t
  use inputs
  use con_data, only: time,dt, dz
  use pars, only: flg_npz, nscl,chem0d,flg_alk, k_ext
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Parameters

  ! LOCATED IN MODULES/PAR.F
  ! flg_reaction : Flag to turn on or off the reaction models

  ! LOCATED IN TRACER/TRACERBC.F90
  ! rmodel       : An array of models to be used in the reactive tracers. Models are
  !                separated out in the react_src function
  ! rdorg        : An array of 0 or 1, 0 = decaying reaction, 1 = growing reaction
  ! rpartner     : An array dictating what tracers are coupled to each other
  ! tau          : The timescale to be used.

  ! flg_debug    : Write a debug file
  integer, parameter :: flg_debug = 0

contains

  ! REACT_SRC: calculate the scalar reaction source term for a given scalar
  !            and point. This is called in rhs_scl for each scalar.
  function react_src(ix,iy,iscl, iz)
    ! ix, iy, iscl, iz (in): location of interest
    real, dimension(nscl-1) :: react_src
    integer, intent(in) :: ix, iy, iscl, iz

    integer :: zzi, i
    real, dimension(0:nscl-2) :: co2, co2tmp
    real :: dtst, temper, K1star, K2star, Kwater, Kboron, Rgas, salt
    real :: h, t_rkc, t_next, t_end, task, k_ext
    integer :: steps

    !!!!!!!!!!!!!
    !! Zeebe et al. 2001 Carbonate Chemistry
    !!!!!!!!!!!!!

    t_rkc  = time
    t_end  = time + dt*0.5
    task   = 1

    ! c(0) = Carbon Dioxide, [CO2], t(ix,iy,2,iz)
    ! c(1) = Bicarbonate, [HCO3-], t(ix,iy,3,iz)
    ! c(2) = Carbonate, [CO32-], t(ix,iy,4,iz)
    ! c(3) = Boric Acid, [B(OH)3], t(ix,iy,5,iz)
    ! c(4) = Tetrahydroxyborate, [B(OH)4-], t(ix,iy,6,iz)
    ! c(5) = Hydrogen Ion, [H+], t(ix,iy,7,iz)
    ! c(6) = Hydroxide, [OH-], t(ix,iy,8,iz)
    
    co2(0:nscl-2)=t(ix, iy, 2:nscl, iz)
    

    if(chem0d == 1) then
      temper = iTsurf
    else
      temper = t(ix,iy,1,iz)
    end if

    co2tmp = intDriver(t_rkc, t_end, co2, temper, iz)
   
    co2=co2tmp

    do i = 0,nscl-2
       react_src(i+1) = co2(i)
    enddo

  end function react_src

  function intDriver(t_rkc, t_end, yGlobal, temper,iz)
    real, intent(in) :: t_rkc
    real, intent(in) :: t_end, temper
    real, intent(in), dimension(0:nscl-2) :: yGlobal
    real, dimension(0:nscl-2) :: intDriver
    real, dimension(0:nscl-2) :: yLocal, yLocal2
    real, dimension(0:4+nscl-1) :: workLocal
    integer i
    integer, intent(in) :: iz

    workLocal(:) = 0.0

    do i = 0,nscl-2
       yLocal(i)    = yGlobal(i)
    enddo

    yLocal2 = rkc_driver(t_rkc, t_end, workLocal, yLocal, temper, iz)


    do i = 0,nscl-2
       yLocal(i)      = yLocal2(i)
       intDriver(i)   = yLocal2(i)
    enddo

  end function intDriver

  function rkc_driver(t_rkc2, t_end, work, yLocal, temper, iz)
    ! Driver function for RKC integrator.
    !
    ! t_tkc    the starting time.
    ! t_end    the desired end time.
    ! task     0 to take a single integration step, 1 to integrate to tEnd.
    ! work     Real work array, size 3.
    ! yLocal   Dependent variable array, integrated values replace initial conditions.

    real, intent(in) :: t_rkc2
    real, intent(in) :: t_end, temper
    real, intent(inout), dimension(0:nscl-2) :: yLocal
    real, dimension(0:nscl-2) :: rkc_driver
    real, intent(inout), dimension(0:4+nscl-1) :: work
    real, dimension(0:nscl-2) :: y_n, F_n, temp_arr, temp_arr2
    integer nstep, m_max, i, m
    real abs_tol, rel_tol, UROUND, hmax, hmin, err, est
    real fac, temp1, temp2, t_rkc
    integer, intent(in) :: iz

    t_rkc = t_rkc2
    nstep   = 0
    rel_tol = 1.0e-6
    abs_tol = 1.0e-10
    UROUND  = 2.22e-16
    m_max   = nint(sqrt(rel_tol / (10.0 * UROUND)))
    hmax    = abs(t_end - t_rkc)
    hmin    = 10.0 * UROUND * max(abs(t_rkc), hmax)

    if(m_max < 2)then
       m_max = 2
    endif

    do i = 0,nscl-2
       y_n(i) = yLocal(i)
    enddo
    

    ! calculate F_n for initial y
    F_n = dydt(t_rkc, y_n, temper, iz)
   
    ! load initial estimate for eigenvector
    if(work(2) < UROUND) then
       do i = 0,nscl-2
          work(4+i) = F_n(i)
       enddo
    endif

    do while (t_rkc < t_end)
       ! use time step stored in work(3)

       ! estimate Jacobian spectral radius
       ! only if 25 steps passed
       ! spec_rad = work(4)
       temp_arr(:)  = 0.0
       temp_arr2(:) = 0.0
       err          = 0.0

       if(mod(nstep,25) == 0)then
          work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4), temp_arr2, temper, iz)
       endif

       ! first step, estimate step size
       if(work(2) < UROUND)then
          work(2) = hmax
          if((work(3) * work(2)) > 1.0)then
             work(2) = 1.0/work(3)
          endif
          work(2) = max(work(2), hmin)

          do i = 0,nscl-2
             temp_arr(i) = y_n(i) + (work(2) * F_n(i))
          enddo
          
          temp_arr2 = dydt(t_rkc + work(2), temp_arr, temper, iz)
          err = 0.0
          do i = 0,nscl-2
             est = (temp_arr2(i) - F_n(i)) / (abs_tol + rel_tol * abs(y_n(i)))
             err = err + est*est
          enddo
          err = work(2) * sqrt(err/real(nscl-2))

          if((0.1 * work(2)) < (hmax * sqrt(err)))then
             work(2) = max((0.1 * work(2)) / sqrt(err), hmin)
          else
             work(2) = hmax
          endif
       endif

       ! check if last step
       if((1.1 * work(2)) .ge. abs(t_end - t_rkc))then
          work(2) = abs(t_end - t_rkc)
       endif

       ! calculate number of steps
       m = 1 + nint(sqrt(1.54 * work(2) * work(3) + 1.0))

       if(m > m_max)then
          m = m_max
          work(2) = real((m*m - 1) / (1.54*work(3)))
       endif

       hmin = 10.0 * UROUND * max(abs(t_rkc), abs(t_rkc + work(2)))

       ! perform tentative time step
       yLocal = rkc_step(t_rkc, work(2), y_n, F_n, m, temper, iz)
       
       ! calculate F_np1 with tenative y_np1
       temp_arr = dydt(t_rkc + work(2), yLocal, temper, iz)
       
       ! estimate error
       err = 0.0
       do i = 0,nscl-2
          est = 0.0
          est = 0.8 * (y_n(i) - yLocal(i)) + 0.4 * work(2) * (F_n(i) + temp_arr(i))
          est = est / (abs_tol + rel_tol * max(abs(yLocal(i)), abs(y_n(i))))
          err = err + est*est
       enddo
       err = sqrt(err / 7.0)

       if (err > 1.0) then
          ! error too large, step is rejected

          ! select smaller step size
          work(2) = 0.8 * work(2) / (err**(1.0/3.0))

          ! reevaluate spectral radius
          work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4), temp_arr2, temper, iz)
          !CALL NPZdebug(y_n(4), y_n(5), y_n(6), iz, 'a rkc_spec_rad2 in rkc_driver')
       else
          ! step accepted
          t_rkc = t_rkc + work(2)
          nstep = nstep + 1

          fac   = 10.0
          temp1 = 0.0
          temp2 = 0.0
          if(work(1) < UROUND)then
             temp2 = err**(1.0/3.0)
             if(0.8 < (fac * temp2))then
                fac = 0.8 /  temp2
             endif
          else
             temp1 = 0.8 * work(2) * (work(0)**(1.0/3.0))
             temp2 = work(1) * (err**(2.0/3.0))
             if(temp1 < (fac * temp2))then
                fac = temp1 / temp2
             endif
          endif

          ! set "old" values to those for current time step
          work(0) = err
          work(1) = work(2)

          do i = 0,nscl-2
             y_n(i) = yLocal(i)
             F_n(i) = temp_arr(i)
          enddo

          ! store next time step
          work(2) = work(2) * max(0.1, fac)
          work(2) = max(hmin, min(hmax, work(2)))

       endif
    enddo

    do i = 0,nscl-2
       rkc_driver(i) = yLocal(i)
    enddo

  end function rkc_driver

  real function rkc_spec_rad(t_rkc, hmax, yLocal, F, v, Fv, temper, iz)
    ! Function to estimate spectral radius.
    !
    ! t_rkc    the time.
    ! hmax     Max time step size.
    ! yLocal   Array of dependent variable.
    ! F        Derivative evaluated at current state
    ! v
    ! Fv

    real, intent(in) :: t_rkc
    real, intent(in) :: hmax, temper
    real, intent(inout), dimension(0:nscl-2) :: v, Fv, F
    real, intent(in), dimension(0:nscl-2) :: yLocal
    integer itmax, i, iter, ind
    real UROUND, small, nrm1, nrm2, dynrm, sigma
    integer, intent(in) :: iz

    UROUND  = 2.22e-16
    itmax   = 50
    small   = 1.0 / hmax
    nrm1    = 0.0
    nrm2    = 0.0
    sigma   = 0.0

    do i = 0,nscl-2
       nrm1 = nrm1 + yLocal(i) * yLocal(i)
       nrm2 = nrm2 + v(i) * v(i)
    enddo
    nrm1 = sqrt(nrm1)
    nrm2 = sqrt(nrm2)

    if((nrm1 .ne. 0.0) .and. (nrm2 .ne. 0.0))then
       dynrm = nrm1 * sqrt(UROUND)
       do i = 0,nscl-2
          v(i) = yLocal(i) + v(i) * (dynrm / nrm2)
       enddo
    elseif(nrm1 .ne. 0.0)then
       dynrm = nrm1 * sqrt(UROUND)
       do i = 0,nscl-2
          v(i) = yLocal(i) * (1.0 + sqrt(UROUND))
       enddo
    elseif(nrm2 .ne. 0.0)then
       dynrm = UROUND
       do i = 0,nscl-2
          v(i) = v(i) * (dynrm / nrm2)
       enddo
    else
       dynrm = UROUND
       do i = 0,nscl-2
          v(i) = UROUND
       enddo
    endif

    ! now iterate using nonlinear power method
    sigma = 0.0
    do iter = 1,itmax
       Fv = dydt(t_rkc, v, temper, iz)

       nrm1 = 0.0
       do i = 0,nscl-2
          nrm1 = nrm1 + ((Fv(i) - F(i)) * (Fv(i) - F(i)))
       enddo
       nrm1  = sqrt(nrm1)
       nrm2  = sigma
       sigma = nrm1 / dynrm
       if((iter .ge. 2) .and. (abs(sigma - nrm2) .le. (max(sigma, small) * 0.01)))then
          do i = 0,nscl-2
             v(i) = v(i) - yLocal(i)
          enddo
          rkc_spec_rad = 1.2 * sigma
       endif

       if(nrm1 .ne. 0.0)then
          do i = 0,nscl-2
             v(i) = yLocal(i) + ((Fv(i) - F(i)) * (dynrm / nrm1))
          enddo
       else
          ind = mod(iter, int(nscl-1))
          v(ind) = yLocal(ind) - (v(ind) - yLocal(ind))
       endif
       !CALL NPZdebug(v(4), v(5), v(6), iz, 'end of rkc_spec_rad')
       
    enddo

    rkc_spec_rad = 1.2 * sigma

  end function rkc_spec_rad

  function rkc_step(t_rkc, h, y_0, F_0, s, temper, iz)
    ! Function to take a single RKC integration step
    !
    ! t_rkc    the starting time.
    ! h        Time-step size.
    ! y_0      Initial conditions.
    ! F_0      Derivative function at initial conditions.
    ! s        number of steps.
    ! rkc_step Integrated variables

    real, intent(in) :: t_rkc
    real, intent(in) :: h, temper
    real, intent(inout), dimension(0:nscl-2) :: y_0, F_0
    integer, intent(in) :: s
    real, dimension(0:nscl-2) :: rkc_step
    real, dimension(0:nscl-2) :: y_j
    real w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
    real c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
    real zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j
    real, dimension(0:nscl-2) :: y_jm1, y_jm2
    integer i, j
    integer, intent(in) :: iz

    w0    = 1.0 + 2.0 / (13.0 * real(s * s))
    temp1 = (w0 * w0) - 1.0
    temp2 = sqrt(temp1)
    arg   = real(s) * log(w0 + temp2)
    w1    = sinh(arg) * temp1 / (cosh(arg) * real(s) * temp2 - w0 * sinh(arg))

    b_jm1 = 1.0 / (4.0 * (w0 * w0))
    b_jm2 = b_jm1

    ! calculate y_1
    mu_t = w1 * b_jm1
    do i = 0,nscl-2
       y_jm2(i) = y_0(i)
       y_jm1(i) = y_0(i) + (mu_t * h * F_0(i))
    enddo

    c_jm2 = 0.0
    c_jm1 = mu_t
    zjm1 = w0
    zjm2 = 1.0
    dzjm1 = 1.0
    dzjm2 = 0.0
    d2zjm1 = 0.0
    d2zjm2 = 0.0

    do j = 2,s

       zj = 2.0 * w0 * zjm1 - zjm2
       dzj = 2.0 * w0 * dzjm1 - dzjm2 + 2.0 * zjm1
       d2zj = 2.0 * w0 * d2zjm1 - d2zjm2 + 4.0 * dzjm1
       b_j = d2zj / (dzj * dzj)
       gamma_t = 1.0 - (zjm1 * b_jm1)

       nu = -b_j / b_jm2
       mu = 2.0 * b_j * w0 / b_jm1
       mu_t = mu * w1 / w0

       ! calculate derivative, use y array for temporary storage
       y_j = dydt(t_rkc + (h * c_jm1), y_jm1, temper, iz)

       do i = 0,nscl-2 
          y_j(i) = (1.0 - mu - nu) * y_0(i) + (mu * y_jm1(i)) + (nu * y_jm2(i)) &
               + h * mu_t * (y_j(i) - (gamma_t * F_0(i)))
       enddo
       
       c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)

       if(j < s)then
          do i = 0,nscl-2
             y_jm2(i) = y_jm1(i)
             y_jm1(i) = y_j(i)
          enddo
       endif

       c_jm2  = c_jm1
       c_jm1  = c_j
       b_jm2  = b_jm1
       b_jm1  = b_j
       zjm2   = zjm1
       zjm1   = zj
       dzjm2  = dzjm1
       dzjm1  = dzj
       d2zjm2 = d2zjm1
       d2zjm1 = d2zj
    enddo

    do i = 0,nscl-2
       rkc_step(i) = y_j(i)
    enddo

  end function rkc_step

  function dydt(t_rkc, y, temper, iz)
   real, intent(in),  dimension(0:nscl-2) :: y
    real, dimension(0:nscl-2) :: dydt, dy
    real, dimension(nscl-1) :: c
    real, intent(in) :: t_rkc, temper
    integer, intent(in) ::  ix, iy, iz
    !     tracers
    real:: P_c, P_n, P_p, P_l, Z_c, Z_n, Z_p, DOM_c, DOM_n, DOM_p, POM_c, POM_n, POM_p, Oxy, Ph, Nit, Am
    !     ratios of tracers
    real:: PnPc, PpPc, PlPc, ZnZc, ZpZc
    !     conversion constants and misc
    real:: E2W = 0.217, CM2M = 0.01, SEC_PER_DAY = 86400.0, HOURS_PER_DAY = 24.0
    real:: p_small = 1.0e-10, pi = 3.14159265359
    !     temperature forcing terms
    real:: Tsum = 25.0,Twin = 8.0, dyear, dfrac
    !     salinity forcing terms
    real:: Salt,Ssum=36.9,Swin=36.9
    !     wind forcing terms
    real:: Wind,Wsum=5.0,Wwin=5.0,C1,C2,C3,C4,O2SCHMIDT
    !     short wave irradiance forcing terms
    real:: Qs,Qsum=30.0,Qwin=30.0,light,declination,dlength,daytime
    !     phytoplankton term coefficients
    real:: Q10p, r_P0, bP, d_P0, h_Pnp, d_Pchi, h_Pchi, betaP, gammaP, h_Pn, n_Pmin, n_Popt
    real:: n_Pmax, p_Pmin, p_Popt, p_Pmax, d_Pchl, alpha_chl0, theta_chl0, c_P
    real:: a_Pph, a_Pn, extra_n, extra_p
    !     zooplankton term coefficients
    real:: Q10z, r_Z0, del_ZP, del_ZZ, bZ, etaZ, betaZ, n_Zopt, p_Zopt, v_Z, d_Zdns, gammaZ
    real:: d_Z0, fZ0, F_c, I_Zc, I_Zn, I_Zp
    !     nutrient and oxygen term coefficients
    real:: Lambda_N4nit, Q10n, h_o, c_POM, eo, abt, h, cxoO2, pschmidt, ScRatio, kun, O2AIRFlux
    !     light attenuation terms
    real,dimension(0:nnz):: EIR, xEPS
    real:: IRR, E_PAR, E_PARk, epsilon_PAR, lambda_w, fEP
    !     temperature effects on biology and chemistry
    real:: fTP, fTZ, fTn
    !     
    real:: epsilon_Pnp, chi_lys, f_Pn, f_Pp, f_Pnp, GP, nu_P, epsilon_P_N, epsilon_PAm, hPnp_fPnp, hPn_N4
    real:: epsilon_Pn, epsilon_Pp, rho_chl, switch_n, switch_p, Gam_Zn, Gam_Zp, QZc, Qzn, Qzp
    real :: N_upt, P_upt, Omega_n, Omega_c, epsilonZc, epsilonZn, epsilonZp, hZF, muZ, v_NZ, v_PZ
    !     parameterization terms for bacteria
    real:: alpha_DOM_O3, alpha_POM_O3, alpha_DOM_N1, alpha_POM_N1, alpha_DOM_N4, alpha_POM_N4
    !     reaction source terms
    real:: dPcdt_gppO, dPcdt_rspO, dPcdt_lysDOMc, dPcdt_lysPOMc, dPcdt_exuDOMc, dPcdt_prdZc
    real:: dPndt_uptN, dPndt_uptAm, dPndt_lysDOMn, dPndt_lysPOMn
    real:: dPpdt_uptPh, dPpdt_lysDOMp, dPpdt_lysPOMp
    real:: dPldt_syn, dPldt_loss
    real:: dZcdt_prdPc, dZcdt_prdZc, dZcdt_rspO, dZcdt_relDOMc, dZcdt_relPOMc
    real:: dZndt_relDOMn, dZndt_relPOMn, dZndt_relAm
    real:: dZpdt_relDOMp, dZpdt_relPOMp, dZpdt_relPh
    real:: dOdt_wind, dAmdt_nitN, dPhdt_denitN2
    integer i, zzi
    do i = 0,nscl-2
       c(i+1) = y(i)
    enddo

    P_c=c(1)
    P_n=c(2)
    P_p=c(3)
    P_l=c(4)
    Z_c=c(5)
    Z_n=c(6)
    Z_p=c(7)
    DOM_c=c(8)
    DOM_n=c(9)
    DOM_p=c(10)
    POM_c=c(11)
    POM_n=c(12)
    POM_p=c(13)
    Oxy=c(14)
    Ph=c(15)
    Nit=c(16)
    Am=c(17)

    ! Set Tracer Ratio Terms
    if(P_n.gt.0.0) PnPc = P_n/P_c
    if(P_c.gt.0.0) PpPc = P_p/P_c
    if(P_l.gt.0.0) PlPc = P_l/P_c
    if(Z_n.gt.0.0) ZnZc = Z_n/Z_c
    if(Z_c.gt.0.0) ZpZc = Z_p/Z_c

    ! Temperary temperature (sin wave) - should be replaced with Temper = t(ix,iy,1,iz)-273.15 (degC)
    dyear         = (time)/SEC_PER_DAY + 0.5
    dfrac         = MOD((time),SEC_PER_DAY)/SEC_PER_DAY
!    Temper        = (Tsum+Twin)/2.0 - ((Tsum-Twin)/2.0)*cos((dyear)*(pi/180.0)) - 0.5*cos(2.0*pi*dfrac)
    ! Salinity
    Salt          = (Ssum+Swin)/2.0 - ((Ssum-Swin)/2.0)*cos((dyear+0.5)*(pi/180.0))
    ! Wind
    Wind          = (Wsum+Wwin)/2.0 - ((Wsum-Wwin)/2.0)*cos((dyear+(dfrac - 0.5))*(pi/180.0))
    ! Short wave irradiance flux - (W/m2)
    light         = (Qsum+Qwin)/2.0 - ((Qsum-Qwin)/2.0)*cos(dyear*(pi/180.0))
    declination   = -0.406*cos(2.0*pi*int(dyear)/360.0)
    dlength       = (acos(-tan(declination)*tan(45.0*(pi/180.0)))/pi*24.0)/2.0
    daytime       = abs(dfrac*24.0-12.0)
    if(daytime.lt.dlength)then
       daytime     = daytime/dlength*pi
       Qs          = light*cos(daytime)+light
    else
       Qs          = 0.0
    end if
    
    ! Phytoplankton coefficients
    Q10p          = 2.00               ! Characteristic Q10 coefficient for phytoplankton
    r_P0          = 1.60/SEC_PER_DAY   ! Maximum specific photosynthetic rate
    bP            = 0.05/SEC_PER_DAY   ! Basal specific respiration rate
    d_P0          = 0.05/SEC_PER_DAY   ! Maximum specific nutrient-stress lysis rate
    h_Pnp         = 0.1                ! Nutrient stress threshold
    betaP         = 0.05               ! Excreted fraction of primary production
    gammaP        = 0.05               ! Activity respiration fraction
    a_Pn          = 0.025/SEC_PER_DAY  ! Specific affinity constant for N
    h_Pn          = 1.0                ! Half saturation constant for ammonium uptake preference over ammonium
    n_Pmin        = 6.87e-3            ! Minimum nitrogen quota
    n_Popt        = 1.26e-2            ! Optimal nitrogen quota
    n_Pmax        = 2.5*n_Popt         ! Maximum nitrogen quota
    a_Pph         = 0.0025/SEC_PER_DAY ! Specific affinity constant for P
    p_Pmin        = 4.29e-4            ! Minimum phosphorus quota
    p_Popt        = 7.86e-4            ! Optimal phosphorus quota
    p_Pmax        = 2.0*p_Popt         ! Maximum phosphorus quota
    alpha_chl0    = 1.52e-5            ! Maximum light utilization coefficient
    theta_chl0    = 0.016              ! Maximum chl:C quotum
    c_P           = 0.03               ! Chl-specific light absorption coefficient
    epsilon_Pn    = n_Pmin/(P_n/P_c)   ! 
    epsilon_Pp    = p_Pmin/(P_p/P_c)   ! 
    ! Mesozooplankton coefficients
    Q10z          = 2.0                ! Characteristic Q10 coefficient for mesozooplankton
    r_Z0          = 2.0/SEC_PER_DAY    ! Potential specific growth rate
    del_ZP        = 1.0                ! Availability of phytoplankton P to zooplankton Z
    del_ZZ        = 0.0                ! Availability of zooplankton Z to zooplankton Z
    bZ            = 0.02/SEC_PER_DAY   ! Basal specific respiration rate
    etaZ          = 0.5                ! Assimilation efficency
    betaZ         = 0.25               ! Excreted fraction of uptake (faeces production)
    epsilonZc     = 0.6                ! Partition between dissolved and particulate excretion of C
    epsilonZn     = 0.72               ! Partition between dissolved and particulate excretion of N
    epsilonZp     = 0.832              ! Partition between dissolved and particulate excretion of P
    hZF           = 200.0              ! Half-saturation food concentration for Type II
    muZ           = 50.0               ! Half-saturation food concentration for preference factor
    n_Zopt        = 0.0126             ! Maximum nitrogen quota
    p_Zopt        = 7.86e-4            ! Maximum phosphorus quota
    d_Z0          = 0.00/SEC_PER_DAY   ! Specific mortality rate
    d_Zdns        = 0.25/SEC_PER_DAY   ! Density-dependent specific mortality rate
    v_PZ          = 1.0/SEC_PER_DAY
    v_NZ          = 1.0/SEC_PER_DAY
    ! Inorganic nutrient coefficients
    Omega_c       = 12.0               ! Unit conversion factor and stoichiometric coefficient
    Omega_n       = 2.0                ! Stoichiometric coefficient for nitrification reaction
    Lambda_N4nit  = 0.01/SEC_PER_DAY   ! Specific nitrification rate at 10 degC
    Q10n          = 2.0                ! Q10 factor for nitrification/denitrification reaction
    h_o           = 10.0               ! Half saturation for chemical processes
    fZ0           = (max(p_small,Oxy)**3)/(max(p_small,Oxy)**3 + h_o**3)  !
    eo            = max(p_small,Oxy)/(max(p_small,Oxy)+h_o)
    ! Wind airation coefficients
    C1            = 1953.4             ! coefficient for oxygen saturation calculation
    C2            = 128.0              ! coefficient for oxygen saturation calculation
    C3            = 3.9918             ! coefficient for oxygen saturation calculation
    C4            = 0.050091           ! coefficient for oxygen saturation calculation
    O2SCHMIDT     = 660.0              ! Schmidt number ratio constant
    ! Irradiance coefficients
    epsilon_PAR   = 0.40               ! Coefficient determining the fraction of PAR in QS
    lambda_w      = 0.0435             ! Background extinction of water.
    ! DOM/POM coefficients 
    c_POM         = 0.1e-3             ! C-specific attenuation coefficient of particulate detritus
    alpha_DOM_O3  = 0.05/SEC_PER_DAY   ! Parameterization of DOM nutrient cycling to nutrient sink
    alpha_POM_O3  = 0.50/SEC_PER_DAY   ! Parameterization of POM nutrient cycling to nutrient sink
    alpha_DOM_N1  = 0.05/SEC_PER_DAY   ! Parameterization of DOM nutrient cycling to inorganic nitrate reservior
    alpha_POM_N1  = 0.50/SEC_PER_DAY   ! Parameterization of POM nutrient cycling to inorganic nitrate reservior
    alpha_DOM_N4  = 0.05/SEC_PER_DAY   ! Parameterization of DOM nutrient cycling to inorganic ammonium reservior
    alpha_POM_N4  = 0.50/SEC_PER_DAY   ! Parameterization of POM nutrient cycling to inorganic ammonium reservior
    
!!! Calculation of Terms !!!
    ! Calculate light attenuation due to water and suspended matter
    E_PARk       = (epsilon_PAR*Qs/E2W)                                                                          ! [eq. 1.1.8] - Photosynthetically available radiation
    EIR(0)    = E_PARk
    xEPS(0)   = lambda_w                                     ! [eq. 1.1.7]
    do zzi = izs,ize
!       if (zzi.gt.1) then ! Calculate vertical extinction of light
          xEPS(zzi) = xEPS(zzi-1) + abs(zl/real(nnz))*((c_POM*t(ix,iy,12,zzi-1) + c_P*t(ix,iy,5,zzi-1) &
               + c_POM*t(ix,iy,12,zzi) + c_P*t(ix,iy,5,zzi))/2.0)
          EIR(zzi)  = (EIR(zzi-1)/(xEPS(zzi)*abs(zl/real(nnz))))*(1.0-exp(-xEPS(zzi)*abs(zl/real(nnz))))
          !       EIR(zzi) = EIR(zzi-1)*exp(-xEPS(zzi)*abs(zl/real(nnz)))
!       end if
    end do
!    if(it.eq.1)then
       if(ix.eq.10)then
          if(iy.eq.10)then
             if(iz.eq.10)then
                WRITE(*,*) 'xEPS = ', xEPS
                WRITE(*,*) 'EIR = ', EIR
                WRITE(*,*) 'izs = ',izs
                WRITE(*,*) 'ize = ',ize
             endif
          endif
       endif
!    endif
    
    E_PAR         = EIR(iz)!/(xEPS(iz)*abs(zl/real(nnz))))*(1.0-exp(-xEPS(iz)*abs(zl/real(nnz))))                                                      ! [eq. 1.1.6]  - Integrated PAR over the level depth
    IRR           = max(p_small,E_PAR*SEC_PER_DAY)                                                                 ! [eq. 1.1.5]  - unit conversion of PAR
    fEP           = 1.0 - exp(-((IRR*alpha_chl0)/(r_P0))*(PlPc))                                                   ! [eq. 1.1.4]  - light regulation factor
    
    ! Calcualte Phytoplankton Carbon Terms
    fTP           = Q10p**((Temper-15.0)/15.0)                                                                         ! [eq. 1.1.3]  - temperature regulation for phytoplankton
    fTP           = max(0.0,fTP)
    fTZ           = Q10z**((Temper-15.0)/15.0)                                                                         ! [eq. 2.1.3]  - temperature regulation for mesozooplankton
    fTZ           = max(0.0,fTZ)
    fTn           = Q10n**((Temper-15.0)/15.0)                                                                         ! [eq. 5.3.3]  - temperature regulation for (de)nitrification processes
    fTn           = max(0.0,fTn)
    f_Pn          = min(1.0,max(p_small,(PnPc - n_Pmin)/(n_Popt - n_Pmin)))                                        ! [eq. 1.1.12] - regulation factor for nitrogen limitation
    f_Pp          = min(1.0,max(p_small,(PpPc - p_Pmin)/(p_Popt - p_Pmin)))                                        ! [eq. 1.1.13] - regulation factor for phosphorous limitation
    f_Pnp         = min(f_Pn,f_Pp)                                                                                 ! [eq. 1.1.11] - regulation factor for multiple nutrient limitation (N and P)
    F_c           = del_ZP*P_c + del_ZZ*Z_c                                                                        ! [eq. 2.1.4]  - availabilbity of prey P for predator Z
    epsilon_Pnp   = min(1.0,p_Pmin/(PpPc + p_small),n_Pmin/(PnPc + p_small))                                       ! [eq. 1.1.16] - portion of phytoplankton lysis going to  DOC
    hPnp_fPnp     = h_Pnp/(h_Pnp+f_Pnp)
    
    dPcdt_gppO    = fTP*fEP*r_P0*P_c                                                                               ! [eq. 1.1.2]  - gross primary production of phytoplankton
    dPcdt_rspO    = gammaP*dPcdt_gppO + fTP*bP*P_c                                                                 ! [eq. 1.1.14] - respiration (carbon loss to inorganic oxygen reservoir)
    dPcdt_lysDOMc = (1.0 - epsilon_Pnp)*hPnp_fPnp*d_P0*P_c                                                         ! [eq. 1.1.15] - lysis (carbon loss to disolved organic matter)
    dPcdt_lysPOMc = epsilon_Pnp*hPnp_fPnp*d_P0*P_c                                                                 ! [eq. 1.1.20] - lysis (carbon loss to particulate organic matter)
    dPcdt_exuDOMc = (betaP + (1-betaP)*(1.0 - f_Pnp))*dPcdt_gppO                                                     ! Exudation
    dPcdt_prdZc   = ((fTZ*r_Z0*(F_c/(F_c+hZF))*Z_c)/(p_small + F_c))*P_c                                ! [eq. 2.1.2]  - predation (carbon loss to mesozooplankton)
    
    ! Calculate Phytoplankton Nitrogen Terms
    GP            = max(0.0,dPcdt_gppO - dPcdt_rspO - dPcdt_lysDOMc - dPcdt_lysPOMc)                               ! [eq. 1.2.4]  - net primary production
    nu_P          = fTP*r_P0                                                                                       ! [eq. 1.2.5]
    hPn_N4        = h_Pn/(p_small+h_Pn+Am)
    epsilon_P_N   = (a_Pn*hPn_N4*Nit*P_c)/(p_small + a_Pn*Am*P_c + a_Pn*hPn_N4*Nit*P_c)                            ! [eq. 1.2.6]  - portion of phytoplankton nitrogen uptake that comes from nitrate reservior
    epsilon_PAm   = (a_Pn*Am*P_c)/(p_small + a_Pn*Am*P_c + a_Pn*hPn_N4*Nit*P_c)                                    ! [eq. 1.2.7]  - portion of phytoplankton nitrogen uptake that comes from ammonium reservior
    N_upt         = min(a_Pn*hPn_N4*Nit*P_c + a_Pn*Am*P_c,n_Popt*GP + nu_P*(n_Pmax*P_c - P_n))                     ! [eq. 1.2.3]  - nitrogen uptake rate
    
    if(N_upt.ge.0.0)then ! - if nitrogen uptake rate is positive, then uptake is divided between coming from the nitrate and ammonium reservoir
       dPndt_uptN  = epsilon_P_N*N_upt
       dPndt_uptAm = epsilon_PAm*N_upt
       extra_n     = 0.0
    else                 ! - if nitrogen uptake is negative, all nitrogen goes to the DON pool
       dPndt_uptN  = 0.0
       dPndt_uptAm = 0.0
       extra_n     = abs(N_upt)
    endif
    dPndt_lysDOMn = (1.0 - epsilon_Pnp)*hPnp_fPnp*d_P0*P_n                                             ! [eq. 1.2.9] - lysis (nitrogen loss to disolved organic matter)
    dPndt_lysPOMn = epsilon_Pnp*hPnp_fPnp*d_P0*P_n                                                      ! [eq. 1.2.10] - lysis (nitrogen loss to particulate organic matter)
    
    ! Phytoplankton Phosphate Terms
    P_upt         = min(a_Pph*Ph*P_c,p_Popt*GP+nu_P*(p_Pmax*P_c - P_p))                                            ! [eq. 1.3.2] - phosphorous uptake rate
    
    if(P_upt.ge.0.0)then ! - if phosphorous uptake rate is positive, then uptake comes from phosphorous reservoir
       dPpdt_uptPh = P_upt
       extra_p     = 0.0
    else                 ! - if phosphorous uptake is negative, all phosphorous goes to the DOP pool
       dPpdt_uptPh = 0.0
       extra_p     = abs(P_upt)
    endif
    dPpdt_lysDOMp = (1.0 - epsilon_Pp)*hPnp_fPnp*d_P0*P_p                                               ! [eq. 1.3.5] - lysis (phosphorous loss to disolved organic matter)
    dPpdt_lysPOMp = epsilon_Pp*hPnp_fPnp*d_P0*P_p                                                       ! [eq. 1.3.6] - lysis (phosphorous loss to particulate organic matter)
    
    ! Phytoplankton Chlorophyll Terms (currently Geider et al. 1997 formulation - PELAGOS)
    rho_chl       = theta_chl0*min(1.0,((fEP*r_P0*P_c)/(alpha_chl0*IRR*(P_l + p_small))))                          ! [eq. 1.4.3] - dynamical chl:C ratio
    dPldt_syn     = rho_chl*dPcdt_gppO*(1.0-gammaP)*(1.0-(1.0-betaP)*(1.0-f_Pnp)-betaP)-hPnp_fPnp*d_P0*P_l
! - PlPc*(dPcdt_lysDOMc + dPcdt_lysPOMc + dPcdt_rspO)                                  ! [eq. 1.4.2] - chlorophyll synthesis
    dPldt_loss    = PlPc*dPcdt_prdZc                                                                               ! [eq. 1.4.4] - chlorophyll loss
    
    ! Mesoszooplankton Carbon Terms
    dZcdt_prdPc   = dPcdt_prdZc                                                                                    ! [eq. 2.1.2] - uptake of carbon (predation of phytoplankton)
    dZcdt_prdZc   = 0.0
                               !
    ! total carbon predated
    I_Zc = dZcdt_prdPc + dZcdt_prdZc
    ! total nitrogen predated
    I_Zn = ZnZc*dZcdt_prdZc + PnPc*dZcdt_prdPc
    ! total phosphate predated
    I_Zp = ZpZc*dZcdt_prdZc + PpPc*dZcdt_prdPc
    
    dZcdt_rspO    = (1.0 - etaZ - betaZ)*I_Zc + bZ*fTZ*Z_c                                                         ! [eq. 2.1.5] - respiration (carbon loss to inorganic oxygen reservoir)
    dZcdt_relDOMc = epsilonZc*(betaZ*I_Zc + (d_Z0 + d_Zdns*(1.0 - fZ0))*fTZ*Z_c)                                     ! [eq. 2.1.6] - excretion/egestion/mortality (carbon loss to disolved organic matter, 0 for mesozooplankton)
    dZcdt_relPOMc = (1.0-epsilonZc)*dZcdt_relDOMc                                                                  ! [eq. 2.1.7] - excretion/egestion/mortality (carbon loss to particulate organic matter)
    
    dZndt_relDOMn = epsilonZn*(betaZ*I_Zn + ZnZc*(d_Z0 + d_Zdns*(1.0 - fZ0))*fTZ*Z_c)                               ! [eq. 2.2.3] - excretion/egestion/mortality (nitrogen loss to disolved organic matter, 0 for mesozooplankton)
    dZndt_relPOMn = (1.0-epsilonZn)*dZndt_relDOMn                                                                  ! [eq. 2.2.4] - excretion/egestion/mortality (nitrogen loss to particulate organic matter)
    dZndt_relAm   = v_NZ*max(0.0,ZnZc-n_Zopt)                                                                      ! [eq. 2.2.5] - excretion/egestion/mortality (nitrogen loss to inorganic ammonium reservior)

    dZpdt_relDOMp = epsilonZp*(betaZ*I_Zp + ZpZc*(d_Z0 + d_Zdns*(1.0 - fZ0))*fTZ*Z_c)                            ! [eq. 2.3.3] - excretion/egestion/mortality (phophorous loss to disolved organic matter, 0 for mesozooplankton)
    dZpdt_relPOMp = (1.0-epsilonZp)*dZpdt_relDOMp                                                                  ! [eq. 2.3.4] - excretion/egestion/mortality (phosphorous loss to particulate organic matter)
    dZpdt_relPh   = v_PZ*max(0.0,ZpZc-p_Zopt)                                                                      ! [eq. 2.3.5] - xcretion/egestion/mortality (phosphorous loss to inorganic phosphate reservior)
    
    ! POM/DOM/Free Nutrient Terms
    dAmdt_nitN    = max(0.0,Lambda_N4nit*fTn*eo*Am)                                                                ! [eq. 5.3.2] - nitrification of ammonium
    
    ! Oxygen airation by wind
    abt           = (Temper+273.15)/100.0
    h             = -173.4292 + (249.6339/abt) + (143.3483*log(abt))-(21.8492*abt)
    h             = h + Salt*(- 0.033096 + 0.014259*abt - 0.0017*(abt**2.0))
    h             = exp(h)
    cxoO2         = h*44.661
    pschmidt      = (C1 - C2*Temper + C3*(Temper**2) - C4*(Temper**3.0))
    ScRatio       = O2SCHMIDT/pschmidt
    kun           = (0.31*(Wind**2.0))*sqrt(ScRatio)*CM2M*HOURS_PER_DAY/SEC_PER_DAY
    O2AIRFlux     = kun*(cxoO2 - Oxy)
    dOdt_wind     = O2AIRFlux/abs(z(iz))
    
    ! This term is then added to the right-hand-side of each tracer's PDE in sr. rhs_scl.f
    ! Phyto Carbon Biochemical Flux [eq. 1.1.1]
    dy(0)  = dPcdt_gppO - dPcdt_rspO - dPcdt_lysDOMc - dPcdt_lysPOMc - dPcdt_exuDOMc - dPcdt_prdZc
    ! Phyto Nitrogen Biochemical Flux [eq. 1.2.1]
    dy(1)  = dPndt_uptN + dPndt_uptAm - dPndt_lysDOMn - dPndt_lysPOMn - PnPc*dPcdt_prdZc  - extra_n
    ! Phyto Phosphate Biochemical Flux [eq. 1.3.1]
    dy(2)  = dPpdt_uptPh - dPpdt_lysDOMp - dPpdt_lysPOMp - PpPc*dPcdt_prdZc  - extra_p
    ! Phyto Chloropyll Biochemical Flux [eq. 1.4.1]
    dy(3)  = dPldt_syn - dPldt_loss
    ! Zoop Carbon Biochemical Flux [eq. 2.1.1]
    dy(4)  = dZcdt_prdPc - dZcdt_rspO - dZcdt_relDOMc - dZcdt_relPOMc
    ! Zoop Nitrogen Biochemical Flux [eq. 2.2.1]
    dy(5)  = PnPc*dZcdt_prdPc - dZndt_relDOMn - dZndt_relPOMn - dZndt_relAm
    ! Zoop Phosphate Biochemical Flux [eq. 2.3.1]
    dy(6)  = PpPc*dZcdt_prdPc - dZpdt_relDOMp - dZpdt_relPOMp - dZpdt_relPh
    ! DOM Carbon Biochemical Flux [eq. 3.1.1]
    dy(7)  = dPcdt_lysDOMc + dPcdt_exuDOMc + dZcdt_relDOMc - alpha_DOM_O3*DOM_c
    ! DOM Nitrogen Biochemical Flux [eq. 3.2.1]
    dy(8)  = dPndt_lysDOMn + dZndt_relPOMn + extra_n - alpha_DOM_N4*DOM_n
    ! DOM Phosphate Biochemical Flux [eq. 3.3.1]
    dy(9)  = dPpdt_lysDOMp + dZpdt_relPOMp + extra_p - alpha_DOM_N1*DOM_p
    ! POM Carbon Biochemical Flux [eq. 4.1.1]
    dy(10)  = dPcdt_lysPOMc + dZcdt_relPOMc - alpha_POM_O3*POM_c
    ! POM Nitrogen Biochemical Flux [eq. 4.2.1]
    dy(11)  = dPndt_lysPOMn + dZndt_relPOMn - alpha_POM_N4*POM_n
    ! POM Phosphate Biochemical Flux [eq. 4.3.1]
    dy(12)  = dPpdt_lysPOMp + dZpdt_relPOMp - alpha_POM_N1*POM_p
    ! Oxygen Biochemical Flux [eq. 5.1.1]
    dy(13)  = dOdt_wind + (dPcdt_gppO - dPcdt_rspO - dZcdt_rspO - alpha_POM_O3*POM_c - alpha_DOM_O3*DOM_c) &
         /Omega_c - dAmdt_nitN*Omega_n
    ! Phosphate Biochemical Flux [eq. 5.2.1]
    dy(14)  = alpha_DOM_N1*DOM_p + alpha_POM_N1*POM_p + dZpdt_relPh - dPpdt_uptPh
    ! Nitrogen Biochemical Flux [eq. 5.3.1]
    dy(15)  = dAmdt_nitN - dPndt_uptN
    ! Ammonium Biochemical Flux [eq. 5.4.1]
    dy(16)  = alpha_DOM_N4*DOM_n + alpha_POM_N4*POM_n + dZndt_relAm - dPndt_uptN - dAmdt_nitN
    
    do i = 0,nscl-2
       dydt(i) = dy(i)
    enddo

  end function dydt

end module reaction