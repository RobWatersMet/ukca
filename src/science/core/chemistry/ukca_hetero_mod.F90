module ukca_hetero_mod
  use yomhook, only : dr_hook, lhook
  use nvtx
  use parkind1, only : jpim, jprb
  implicit none
  real, parameter, private :: rho1 = 1.35
  real, parameter, private :: rho2 = 0.928
  real, parameter, private :: rad1 = 1.0e-4
  real, parameter, private :: rad2 = 10.0e-4
  real, parameter, private :: radsa = 1.0e-5
  character(len = *), parameter, private :: ModuleName = 'UKCA_HETERO_MOD'
  public

  contains
  subroutine ukca_hetero(n_points, have_nat, stratflag)
    use asad_mod, only : cdt, f, fpsc1, jpcspf, jphk, nhrkx, p, peps, rk, sh2o, shno3, specf, sph, sph2o, sphno3, t, tnd, wp, za
    use ukca_config_constants_mod, only : avogadro, boltzmann
    use ukca_constants, only : m_brono2, m_clono2, m_h2o, m_hno3, m_hobr, m_hocl, m_n2o5, pi
    use ukca_config_specification_mod, only : ukca_config
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    logical, dimension(n_points), intent(in) :: have_nat
    logical, dimension(n_points), intent(in) :: stratflag
    logical, save :: gpsa
    logical, save :: gphocl
    logical, save :: gppsc
    logical, save :: gpsimp
    logical :: l_ukca_sulphur
    logical :: l_ukca_presaer
    integer :: js
    integer :: jh
    integer :: jl
    integer :: n_hk
    integer, save :: ih2o = 0
    integer, save :: ihno3 = 0
    integer, save :: ihcl = 0
    integer, save :: iclono2 = 0
    integer, save :: ihocl = 0
    integer, save :: in2o5 = 0
    integer, save :: ihbr = 0
    integer, save :: ibrono2 = 0
    integer, save :: ihobr = 0
    integer, save :: n_clono2_hcl = 0
    integer, save :: n_clono2_h2o = 0
    integer, save :: n_n2o5_h2o = 0
    integer, save :: n_n2o5_hcl = 0
    integer, save :: n_hocl_hcl = 0
    integer, save :: n_brono2_hcl = 0
    integer, save :: n_brono2_h2o = 0
    integer, save :: n_hobr_hcl = 0
    integer, save :: n_hobr_hbr = 0
    integer, save :: n_hocl_hbr = 0
    integer, save :: n_clono2_hbr = 0
    integer, save :: n_brono2_hbr = 0
    integer, save :: n_n2o5_hbr = 0
    logical, save :: first = .true.
    logical, save :: first_pass = .true.
    real, dimension(n_points) :: zp
    real, dimension(n_points) :: zt
    real, dimension(n_points) :: zhno3
    real, dimension(n_points) :: zh2o
    real, dimension(n_points) :: zhcl
    real, dimension(n_points) :: zclono2
    real, dimension(n_points) :: zn2o5
    real, dimension(n_points) :: zhocl
    real, dimension(n_points) :: zbrono2
    real, dimension(n_points) :: zhobr
    real, dimension(n_points) :: zhbr
    real, dimension(n_points) :: psc1
    real, dimension(n_points) :: psc2
    real, dimension(n_points) :: psc3
    real, dimension(n_points) :: psc4
    real, dimension(n_points) :: psc5
    real, dimension(n_points,5) :: hk
    real :: zrate
    real :: zdhcl
    real :: zdhbr
    real :: zfact
    real, dimension(n_points) :: c_cf
    real, dimension(n_points) :: psc1sa
    real, dimension(n_points) :: psc2sa
    real, dimension(n_points,13) :: kpsc
    real, dimension(13) :: gam1
    real, dimension(13) :: gam2
    real, dimension(13) :: gam3
    real, dimension(n_points,7) :: gam_calc
    real, dimension(n_points,8) :: gam3calc
    real, dimension(13) :: mm_arr
    real :: amu
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_HETERO'
    integer :: widx1
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19
    integer :: idx_20
    integer :: idx_21
    integer :: idx_22
    integer :: idx_23
    integer :: idx_24
    integer :: idx_25
    integer :: idx_26
    integer :: idx_27
    integer :: idx_28
    integer :: idx_29
    integer :: idx_30
    integer :: idx_31
    integer :: idx_32
    integer :: idx_33
    integer :: idx_34
    integer :: idx_35
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    gam1(:) = [0.2, 0.004, 0.1, 0.0004, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 7.0]
    gam2(:) = [0.3, 0.3, 0.2, 0.02, 0.03, 0.25, 0.3, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0]
    gam3(:) = [0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mm_arr(:) = [m_clono2, m_clono2, m_hocl, m_n2o5, m_n2o5, m_hobr, m_brono2, m_brono2, m_hobr, m_hocl, m_clono2, m_brono2, m_n2o5]
    if (first_pass) then
      if (first) then
        !NO!$acc kernels
        !NO!$acc loop seq
        do js = 1, jpcspf, 1
          if (specf(js) == 'H2O       ' .OR. specf(js) == 'H2OS      ') then
            ih2o = js
          else
            if (specf(js) == 'HONO2     ') then
              ihno3 = js
            else
              if (specf(js) == 'HCl       ') then
                ihcl = js
              else
                if (specf(js) == 'ClONO2    ') then
                  iclono2 = js
                else
                  if (specf(js) == 'N2O5      ') then
                    in2o5 = js
                  else
                    if (specf(js) == 'HOCl      ') then
                      ihocl = js
                    else
                      if (specf(js) == 'HBr       ') then
                        ihbr = js
                      else
                        if (specf(js) == 'BrONO2    ') then
                          ibrono2 = js
                        else
                          if (specf(js) == 'HOBr      ') then
                            ihobr = js
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        enddo
        !NO!$acc loop gang vector independent
        do jh = 1, jphk, 1
          if (sph(jh,1) == 'H2O       ' .OR. sph(jh,1) == 'H2OS      ') then
            if (sph(jh,2) == 'ClONO2    ') then
              n_clono2_h2o = nhrkx(jh)
            end if
            if (sph(jh,2) == 'N2O5      ') then
              n_n2o5_h2o = nhrkx(jh)
            end if
            if (sph(jh,2) == 'BrONO2    ') then
              n_brono2_h2o = nhrkx(jh)
            end if
          else
            if (sph(jh,1) == 'ClONO2    ') then
              if (sph(jh,2) == 'HCl       ') then
                n_clono2_hcl = nhrkx(jh)
              end if
              if (sph(jh,2) == 'H2O       ') then
                n_clono2_h2o = nhrkx(jh)
              end if
              if (sph(jh,2) == 'H2OS      ') then
                n_clono2_h2o = nhrkx(jh)
              end if
              if (sph(jh,2) == 'HBr       ') then
                n_clono2_hbr = nhrkx(jh)
              end if
            else
              if (sph(jh,1) == 'HCl       ') then
                if (sph(jh,2) == 'ClONO2    ') then
                  n_clono2_hcl = nhrkx(jh)
                end if
                if (sph(jh,2) == 'N2O5      ') then
                  n_n2o5_hcl = nhrkx(jh)
                end if
                if (sph(jh,2) == 'HOCl      ') then
                  n_hocl_hcl = nhrkx(jh)
                end if
                if (sph(jh,2) == 'HOBr      ') then
                  n_hobr_hcl = nhrkx(jh)
                end if
                if (sph(jh,2) == 'BrONO2    ') then
                  n_brono2_hcl = nhrkx(jh)
                end if
              else
                if (sph(jh,1) == 'BrONO2    ') then
                  if (sph(jh,2) == 'HCl       ') then
                    n_brono2_hcl = nhrkx(jh)
                  end if
                  if (sph(jh,2) == 'H2O       ') then
                    n_brono2_h2o = nhrkx(jh)
                  end if
                  if (sph(jh,2) == 'H2OS      ') then
                    n_brono2_h2o = nhrkx(jh)
                  end if
                  if (sph(jh,2) == 'HBr       ') then
                    n_brono2_hbr = nhrkx(jh)
                  end if
                else
                  if (sph(jh,1) == 'HOBr      ') then
                    if (sph(jh,2) == 'HCl       ') then
                      n_hobr_hcl = nhrkx(jh)
                    end if
                    if (sph(jh,2) == 'HBr       ') then
                      n_hobr_hbr = nhrkx(jh)
                    end if
                  else
                    if (sph(jh,1) == 'N2O5      ') then
                      if (sph(jh,2) == 'H2O       ') then
                        n_n2o5_h2o = nhrkx(jh)
                      end if
                      if (sph(jh,2) == 'H2OS      ') then
                        n_n2o5_h2o = nhrkx(jh)
                      end if
                      if (sph(jh,2) == 'HCl       ') then
                        n_n2o5_hcl = nhrkx(jh)
                      end if
                      if (sph(jh,2) == 'HBr       ') then
                        n_n2o5_hbr = nhrkx(jh)
                      end if
                    else
                      if (sph(jh,1) == 'HOCl      ') then
                        if (sph(jh,2) == 'HCl       ') then
                          n_hocl_hcl = nhrkx(jh)
                        end if
                        if (sph(jh,2) == 'HBr       ') then
                          n_hocl_hbr = nhrkx(jh)
                        end if
                      else
                        if (sph(jh,1) == 'HBr       ') then
                          if (sph(jh,2) == 'ClONO2    ') then
                            n_clono2_hbr = nhrkx(jh)
                          end if
                          if (sph(jh,2) == 'BrONO2    ') then
                            n_brono2_hbr = nhrkx(jh)
                          end if
                          if (sph(jh,2) == 'HOBr      ') then
                            n_hobr_hbr = nhrkx(jh)
                          end if
                          if (sph(jh,2) == 'HOCl      ') then
                            n_hocl_hbr = nhrkx(jh)
                          end if
                          if (sph(jh,2) == 'N2O5      ') then
                            n_n2o5_hbr = nhrkx(jh)
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        enddo
        first = .false.
        l_ukca_sulphur = .false.
        l_ukca_presaer = .true.
        gpsa = l_ukca_sulphur .OR. l_ukca_presaer
        gphocl = .true.
        gppsc = .true.
        gpsimp = .false.
        !NO!$acc end kernels
      end if
      first_pass = .false.
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx = 1, n_points, 1
      zp(idx) = p(idx) / 100.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = 1, n_points, 1
      zt(idx_1) = t(idx_1)
    enddo
    !NO!$acc end kernels
    if (ihno3 > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_2 = 1, n_points, 1
        zhno3(idx_2) = f(idx_2 + (LBOUND(f, dim=1) - 1),ihno3)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_3 = 1, n_points, 1
        zhno3(idx_3) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (ih2o > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_4 = 1, n_points, 1
        zh2o(idx_4) = f(idx_4 + (LBOUND(f, dim=1) - 1),ih2o)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_5 = 1, n_points, 1
        zh2o(idx_5) = wp(idx_5 + (LBOUND(wp, dim=1) - 1))
      enddo
      !NO!$acc end kernels
    end if
    if (ihcl > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_6 = 1, n_points, 1
        zhcl(idx_6) = f(idx_6 + (LBOUND(f, dim=1) - 1),ihcl)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_7 = 1, n_points, 1
        zhcl(idx_7) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (iclono2 > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_8 = 1, n_points, 1
        zclono2(idx_8) = f(idx_8 + (LBOUND(f, dim=1) - 1),iclono2)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_9 = 1, n_points, 1
        zclono2(idx_9) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (in2o5 > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_10 = 1, n_points, 1
        zn2o5(idx_10) = f(idx_10 + (LBOUND(f, dim=1) - 1),in2o5)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_11 = 1, n_points, 1
        zn2o5(idx_11) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (ihocl > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_12 = 1, n_points, 1
        zhocl(idx_12) = f(idx_12 + (LBOUND(f, dim=1) - 1),ihocl)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_13 = 1, n_points, 1
        zhocl(idx_13) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (ihobr > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_14 = 1, n_points, 1
        zhobr(idx_14) = f(idx_14 + (LBOUND(f, dim=1) - 1),ihobr)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_15 = 1, n_points, 1
        zhobr(idx_15) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (ihbr > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_16 = 1, n_points, 1
        zhbr(idx_16) = f(idx_16 + (LBOUND(f, dim=1) - 1),ihbr)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_17 = 1, n_points, 1
        zhbr(idx_17) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    if (ibrono2 > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_18 = 1, n_points, 1
        zbrono2(idx_18) = f(idx_18 + (LBOUND(f, dim=1) - 1),ibrono2)
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_19 = 1, n_points, 1
        zbrono2(idx_19) = 0.0
      enddo
      !NO!$acc end kernels
    end if
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (.NOT. (stratflag)) sph2o = 0.0
    call ukca_pscpres(zt(:), zp(:), tnd(1:n_points), zh2o(:), zhno3(:), 1, n_points, n_points, have_nat(:), sph2o(1:n_points))
    if (ihno3 > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_20 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_20,ihno3) = zhno3(idx_20 + (1 - LBOUND(f, dim=1)))
      enddo
      !NO!$acc end kernels
    end if
    if (ih2o > 0) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do idx_21 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_21,ih2o) = zh2o(idx_21 + (1 - LBOUND(f, dim=1)))
      enddo
      !NO!$acc end kernels
    end if
    if (ukca_config%i_ukca_hetconfig == 0) then
      call ukca_calckpsc(za(1:n_points), zt(:), zh2o(:), zhcl(:), zclono2(:), zn2o5(:), zhocl(:), psc1(:), psc2(:), psc3(:), psc4(:), psc5(:), gpsa, gphocl, gppsc, gpsimp, n_points, 1, n_points, cdt)
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (zhcl > peps)
        hk(:, 2) = psc1 / zhcl
        hk(:, 3) = psc5 / zhcl
        hk(:, 5) = psc4 / zhcl
      ELSEWHERE
        hk(:, 2) = 0.0
        hk(:, 3) = 0.0
        hk(:, 5) = 0.0
      END WHERE
      WHERE (zh2o > peps)
        hk(:, 1) = psc2 / zh2o
        hk(:, 4) = psc3 / zh2o
      ELSEWHERE
        hk(:, 1) = 0.0
        hk(:, 4) = 0.0
      END WHERE
      if (n_clono2_h2o > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_22 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_22,n_clono2_h2o) = hk(idx_22 + (LBOUND(hk, dim=1) - LBOUND(rk, dim=1)),1)
        enddo
        !NO!$acc end kernels
      end if
      if (n_clono2_hcl > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_23 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_23,n_clono2_hcl) = hk(idx_23 + (LBOUND(hk, dim=1) - LBOUND(rk, dim=1)),2)
        enddo
        !NO!$acc end kernels
      end if
      if (n_hocl_hcl > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_24 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_24,n_hocl_hcl) = hk(idx_24 + (LBOUND(hk, dim=1) - LBOUND(rk, dim=1)),3)
        enddo
        !NO!$acc end kernels
      end if
      if (n_n2o5_h2o > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_25 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_25,n_n2o5_h2o) = hk(idx_25 + (LBOUND(hk, dim=1) - LBOUND(rk, dim=1)),4)
        enddo
        !NO!$acc end kernels
      end if
      if (n_n2o5_hcl > 0) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_26 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
          rk(idx_26,n_n2o5_hcl) = hk(idx_26 + (LBOUND(hk, dim=1) - LBOUND(rk, dim=1)),5)
        enddo
        !NO!$acc end kernels
      end if
    else
      if (ukca_config%i_ukca_hetconfig == 2) then
        call ukca_br_gamma(1, n_points, n_points, zt(:), zhbr(:), gam_calc)
      end if
      if (ukca_config%i_ukca_hetconfig > 0) then
        call ukca_shi_liquid_aerosol(1, n_points, n_points, zt(:), zh2o(:), zhcl(:), zclono2(:), gam3calc(:,:))
      end if
      if (gppsc) then
        amu = 1.0 / (avogadro * 1000.0)
        if (ukca_config%i_ukca_hetconfig == 2) then
          n_hk = 13
        else
          n_hk = 5
        end if
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(2)
        do idx_27 = LBOUND(kpsc, dim=2), UBOUND(kpsc, dim=2), 1
          do idx_28 = 1, n_points, 1
            kpsc(idx_28,idx_27) = 0.0
          enddo
        enddo
        !NO!$acc loop seq
        do jh = 1, n_hk, 1
          !NO!$acc loop gang vector independent
          do idx_29 = 1, n_points, 1
            c_cf(idx_29) = 0.25 * SQRT((8.0 * boltzmann) / (pi * mm_arr(jh) * amu) * zt(idx_29))
          enddo
          if (gppsc) then
            if (gpsimp) then
              !NO!$acc loop gang vector independent
              do idx_30 = 1, n_points, 1
                kpsc(idx_30,jh) = 4.6e-5 * fpsc1(idx_30)
              enddo
            else
              !NO!$acc loop gang vector independent
              do idx_31 = 1, n_points, 1
                psc1sa(idx_31) = (m_hno3 + 3.0 * m_h2o) * amu * 3.0e5 * shno3(idx_31) / (rho1 * rad1)
              enddo
              !NO!$acc loop gang vector independent
              do idx_32 = 1, n_points, 1
                psc2sa(idx_32) = m_h2o * amu * 3.0e5 * sh2o(idx_32) / (rho2 * rad2)
              enddo
              if (gam1(jh) >= 1.0) then
                !NO!$acc loop gang vector independent
                do idx_33 = 1, n_points, 1
                  kpsc(idx_33,jh) = c_cf(idx_33) * psc1sa(idx_33) * gam_calc(idx_33,INT(gam1(jh)))
                enddo
              else
                !NO!$acc loop gang vector independent
                do idx_34 = 1, n_points, 1
                  kpsc(idx_34,jh) = c_cf(idx_34) * psc1sa(idx_34) * gam1(jh)
                enddo
              end if
              if (gam2(jh) >= 1.0) then
                !NO!$acc loop gang vector independent
                do idx_35 = 1, n_points, 1
                  kpsc(idx_35,jh) = kpsc(idx_35,jh) + c_cf(idx_35) * psc2sa(idx_35) * gam_calc(idx_35,INT(gam2(jh)))
                enddo
              else
                !NO!$acc loop gang vector independent
                do idx_36 = 1, n_points, 1
                  kpsc(idx_36,jh) = kpsc(idx_36,jh) + c_cf(idx_36) * psc2sa(idx_36) * gam2(jh)
                enddo
              end if
            end if
          end if
          if (ukca_config%i_ukca_hetconfig > 0) then
            if (jh <= 3 .OR. jh == 8) then
              !NO!$acc loop gang vector independent
              do idx_37 = 1, n_points, 1
                kpsc(idx_37,jh) = kpsc(idx_37,jh) + c_cf(idx_37) * 100.0 * za(idx_37) * gam3calc(idx_37,jh)
              enddo
            else
              !NO!$acc loop gang vector independent
              do idx_38 = 1, n_points, 1
                kpsc(idx_38,jh) = kpsc(idx_38,jh) + c_cf(idx_38) * 100.0 * za(idx_38) * gam3(jh)
              enddo
            end if
          end if
        enddo
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          zrate = MAX(1.0, cdt * (kpsc(jl,1) * zclono2(jl) + kpsc(jl,5) * zn2o5(jl) + kpsc(jl,3) * zhocl(jl) + kpsc(jl,6) * zhobr(jl) + kpsc(jl,7) * zbrono2(jl)))
          zdhcl = MIN(zrate, zhcl(jl))
          zfact = zdhcl / zrate
          kpsc(jl,1) = zfact * kpsc(jl,1)
          kpsc(jl,5) = zfact * kpsc(jl,5)
          kpsc(jl,3) = zfact * kpsc(jl,3)
          kpsc(jl,6) = zfact * kpsc(jl,6)
          kpsc(jl,7) = zfact * kpsc(jl,7)
        enddo
        !NO!$acc end kernels
        if (ukca_config%i_ukca_hetconfig == 2) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do jl = 1, n_points, 1
            zrate = MAX(1.0, cdt * (kpsc(jl,9) * zhobr(jl) + kpsc(jl,10) * zhocl(jl) + kpsc(jl,11) * zclono2(jl) + kpsc(jl,12) * zbrono2(jl) + kpsc(jl,13) * zn2o5(jl)))
            zdhbr = MIN(zrate, zhbr(jl))
            zfact = zdhbr / zrate
            kpsc(jl,9) = zfact * kpsc(jl,9)
            kpsc(jl,10) = zfact * kpsc(jl,10)
            kpsc(jl,11) = zfact * kpsc(jl,11)
            kpsc(jl,12) = zfact * kpsc(jl,12)
            kpsc(jl,13) = zfact * kpsc(jl,13)
          enddo
          !NO!$acc end kernels
        end if
      end if
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
      WHERE (zhcl > peps)
        rk(:, n_clono2_hcl) = kpsc(:, 1) / zhcl
        rk(:, n_hocl_hcl) = kpsc(:, 3) / zhcl
        rk(:, n_n2o5_hcl) = kpsc(:, 5) / zhcl
        rk(:, n_hobr_hcl) = kpsc(:, 6) / zhcl
        rk(:, n_brono2_hcl) = kpsc(:, 7) / zhcl
      ELSEWHERE
        rk(:, n_clono2_hcl) = 0.0
        rk(:, n_hocl_hcl) = 0.0
        rk(:, n_n2o5_hcl) = 0.0
        rk(:, n_hobr_hcl) = 0.0
        rk(:, n_brono2_hcl) = 0.0
      END WHERE
      WHERE (zh2o > peps)
        rk(:, n_clono2_h2o) = kpsc(:, 2) / zh2o
        rk(:, n_n2o5_h2o) = kpsc(:, 4) / zh2o
        rk(:, n_brono2_h2o) = kpsc(:, 8) / zh2o
      ELSEWHERE
        rk(:, n_clono2_h2o) = 0.0
        rk(:, n_n2o5_h2o) = 0.0
        rk(:, n_brono2_h2o) = 0.0
      END WHERE
      WHERE (zhbr > peps)
        rk(:, n_hobr_hbr) = kpsc(:, 9) / zhbr
        rk(:, n_hocl_hbr) = kpsc(:, 10) / zhbr
        rk(:, n_clono2_hbr) = kpsc(:, 11) / zhbr
        rk(:, n_brono2_hbr) = kpsc(:, 12) / zhbr
        rk(:, n_n2o5_hbr) = kpsc(:, 13) / zhbr
      ELSEWHERE
        rk(:, n_hobr_hbr) = 0.0
        rk(:, n_hocl_hbr) = 0.0
        rk(:, n_clono2_hbr) = 0.0
        rk(:, n_brono2_hbr) = 0.0
        rk(:, n_n2o5_hbr) = 0.0
      END WHERE
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx_39 = LBOUND(sphno3, dim=1), UBOUND(sphno3, dim=1), 1
      sphno3(idx_39) = shno3(idx_39 + (LBOUND(shno3, dim=1) - LBOUND(sphno3, dim=1)))
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_hetero
  subroutine ukca_br_gamma(kstart, kend, kchmlev, t, thbr, gam_br)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kchmlev
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    real, dimension(kchmlev), intent(in) :: t
    real, dimension(kchmlev), intent(in) :: thbr
    real, dimension(kchmlev,7), intent(out) :: gam_br
    real :: tt
    real :: theta_hbr
    integer :: jl
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_BR_GAMMA'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do jl = kstart, kend, 1
      tt = MAX(t(jl), 180.0)
      theta_hbr = 4.14e-10 * thbr(jl) ** 0.88
      gam_br(jl,1) = 5.3e-4 * EXP(1100.0 / tt)
      gam_br(jl,2) = 4.8e-4 * EXP(1240.0 / tt)
      gam_br(jl,3) = theta_hbr / (0.3 * theta_hbr + 2.7)
      gam_br(jl,4) = 0.56 * theta_hbr
      gam_br(jl,5) = 6.6e-3 * EXP(700.0 / tt)
      gam_br(jl,6) = 0.56 * theta_hbr
      gam_br(jl,7) = 0.02 * theta_hbr
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_br_gamma
  subroutine ukca_shi_liquid_aerosol(kstart, kend, kchmlev, t, th2o, thcl, tcnit, gam3arr)
    use ukca_config_constants_mod, only : boltzmann
    use asad_mod, only : peps
    real, parameter :: r_vol = 0.082
    real, parameter :: min_value = 1.0e-15
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kchmlev
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    real, dimension(kchmlev), intent(in) :: t
    real, dimension(kchmlev), intent(in) :: th2o
    real, dimension(kchmlev), intent(in) :: thcl
    real, dimension(kchmlev), intent(in) :: tcnit
    real, dimension(kchmlev,8), intent(out) :: gam3arr
    real :: tt
    real :: p0h2o
    real :: ph2o
    real :: phcl
    real :: pclono2
    real :: aw
    real :: ml_h2so4
    real :: wt_h2so4
    real :: mr_h2so4
    real :: rho_h2so4
    real :: vis_h2so4
    real :: a1
    real :: b1
    real :: c1
    real :: d1
    real :: a2
    real :: b2
    real :: c2
    real :: d2
    real :: z1
    real :: z2
    real :: z3
    real :: mf_h2so4
    real :: a
    real :: t0
    real :: ah
    real :: y1
    real :: y2
    real :: k_h2o
    real :: k_hydr
    real :: k_hcl
    real :: k_hocl
    real :: k_h
    real :: c_clono2
    real :: c_hocl
    real :: h_clono2
    real :: h_hocl
    real :: h_hcl
    real :: s_clono2
    real :: s_hocl
    real :: d_clono2
    real :: d_hocl
    real :: g_h2o
    real :: g_hcl
    real :: g_s
    real :: g_s2
    real :: g_hcl2
    real :: g_clono2
    real :: g_b
    real :: g_s2_b
    real :: gam_hocl_f_hcl
    real :: m_hcl
    real :: l_clono2
    real :: l_hocl
    real :: f_clono2
    real :: f_hocl
    real :: f_hcl
    real :: gam_hocl
    real :: gam_clono2
    real :: small_value
    integer :: jl
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SHI_LIQUID_AEROSOL'
    integer :: idx
    integer :: idx_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do idx = LBOUND(gam3arr, dim=2), UBOUND(gam3arr, dim=2), 1
      do idx_1 = LBOUND(gam3arr, dim=1), UBOUND(gam3arr, dim=1), 1
        gam3arr(idx_1,idx) = 0.0
      enddo
    enddo
    small_value = SQRT(peps)
    !NO!$acc loop gang vector independent
    do jl = kstart, kend, 1
      tt = MAX(t(jl), 185.0)
      p0h2o = EXP(18.452406985 - 3505.1578807 / tt - 330918.55082 / tt ** 2 + 12725068.262 / tt ** 3)
      ph2o = 1.0e4 * th2o(jl) * boltzmann * tt
      phcl = 1.0e4 * thcl(jl) * boltzmann * tt / 1013.0
      pclono2 = 1.0e4 * tcnit(jl) * boltzmann * tt / 1013.0
      if (ABS(p0h2o) < small_value) then
        aw = 0.0
      else
        aw = ph2o / p0h2o
        if (aw < min_value) then
          aw = min_value
        end if
      end if
      if (aw < 0.85) then
        if (aw <= 0.05) then
          a1 = 12.37208932
          b1 = -0.16125516114
          c1 = -30.490657554
          d1 = -2.1133114241
          a2 = 13.455394705
          b2 = -0.1921312255
          c2 = -34.285174607
          d2 = -1.7620073078
        else
          a1 = 11.820654354
          b1 = -0.20786404244
          c1 = -4.807306373
          d1 = -5.1727540348
          a2 = 12.891938068
          b2 = -0.23233847708
          c2 = -6.4261237757
          d2 = -4.9005471319
        end if
      else
        a1 = -180.06541028
        b1 = -0.38601102592
        c1 = -93.317846778
        d1 = 273.88132245
        a2 = -176.95814097
        b2 = -0.36257048154
        c2 = -90.469744201
        d2 = 267.45509988
      end if
      y1 = a1 * aw ** b1 + c1 * aw + d1
      y2 = a2 * aw ** b2 + c2 * aw + d2
      ml_h2so4 = y1 + ((tt - 190.0) * (y2 - y1)) / 70.0
      ml_h2so4 = MAX(ml_h2so4, small_value)
      wt_h2so4 = (9800.0 * ml_h2so4) / (98.0 * ml_h2so4 + 1000.0)
      z1 = 0.12364 - 5.6e-7 * tt ** 2.0
      z2 = -0.02954 + 1.814e-7 * tt ** 2.0
      z3 = 2.343e-3 - 1.487e-6 * tt - 1.324e-8 * tt ** 2.0
      rho_h2so4 = 1.0 + z1 * ml_h2so4 + z2 * ml_h2so4 ** 1.5 + z3 * ml_h2so4 ** 2.0
      mr_h2so4 = (rho_h2so4 * wt_h2so4) / 9.8
      mf_h2so4 = wt_h2so4 / (wt_h2so4 + (100.0 - wt_h2so4) * (98.0 / 18.0))
      if (wt_h2so4 <= 80.0) then
        a = 169.5 + 5.18 * wt_h2so4 - 0.0825 * wt_h2so4 ** 2.0 + 3.27e-3 * wt_h2so4 ** 3.0
        t0 = 144.11 + 0.166 * wt_h2so4 - 0.015 * wt_h2so4 ** 2.0 + 2.18e-4 * wt_h2so4 ** 3.0
        vis_h2so4 = a * tt ** (-1.43) * EXP(448.0 / (tt - t0))
      else
        vis_h2so4 = 1730.14 * tt ** (-1.43) * EXP(448.0 / (tt - 173.006))
      end if
      ah = EXP(60.51 - 0.095 * wt_h2so4 + 0.0077 * wt_h2so4 ** 2.0 - 1.61e-5 * wt_h2so4 ** 3.0 - (1.76 + 2.52e-4 * wt_h2so4 ** 2.0) * tt ** 0.5 + (-805.89 + 253.05 * wt_h2so4 ** 0.076) / tt ** 0.5)
      k_h2o = 1.95e10 * EXP(-2800.0 / tt)
      k_h = 1.22e12 * EXP(-6200.0 / tt)
      k_hydr = k_h2o * aw + k_h * ah * aw
      c_clono2 = 1474.0 * tt ** 0.5
      s_clono2 = 0.306 + 24.0 / tt
      h_clono2 = 1.6e-6 * EXP(4710.0 / tt - s_clono2 * mr_h2so4)
      d_clono2 = 5e-8 * tt / vis_h2so4
      g_h2o = (4 * h_clono2 * r_vol * tt / c_clono2) * (d_clono2 * k_hydr) ** 0.5
      h_hcl = (0.094 - 0.61 * mf_h2so4 + 1.2 * mf_h2so4 ** 2) * EXP(-8.68 + (8515.0 - 10718.0 * mf_h2so4 ** 0.7) / tt)
      m_hcl = h_hcl * phcl
      k_hcl = 7.9e11 * ah * d_clono2 * m_hcl
      l_clono2 = (d_clono2 / (k_hydr + k_hcl)) ** 0.5
      f_clono2 = 1 / TANH(radsa / l_clono2) - l_clono2 / radsa
      g_clono2 = f_clono2 * g_h2o * (1 + k_hcl / k_hydr) ** 0.5
      g_hcl = (g_clono2 * k_hcl) / (k_hcl + k_hydr)
      g_s = 66.12 * EXP(-1374.0 / tt) * h_clono2 * m_hcl
      f_hcl = 1 / (1 + 0.612 * (g_s + g_hcl) * (pclono2 / phcl))
      g_s2 = f_hcl * g_s
      g_hcl2 = f_hcl * g_hcl
      g_b = g_hcl2 + (g_clono2 * k_hydr) / (k_hcl + k_hydr)
      if (ABS(g_s2 + g_b) < small_value) then
        gam3arr(jl,1) = 0.0
        gam3arr(jl,2) = 0.0
      else
        g_s2_b = 1.0 + 1.0 / (g_s2 + g_b)
        if (ABS(g_s2_b) < small_value) then
          gam3arr(jl,1) = 0.0
          gam3arr(jl,2) = 0.0
        else
          gam_clono2 = 1.0 / g_s2_b
          gam3arr(jl,1) = gam_clono2 * ((g_s2 + g_hcl2) / (g_s2 + g_b))
          gam3arr(jl,2) = gam_clono2 - gam3arr(jl,1)
          if (gam3arr(jl,1) < min_value) then
            gam3arr(jl,1) = min_value
          end if
          if (gam3arr(jl,2) < min_value) then
            gam3arr(jl,2) = min_value
          end if
        end if
      end if
      s_hocl = 0.0776 + 59.18 / tt
      h_hocl = 1.91e-6 * EXP(5862.4 / tt - s_hocl * mr_h2so4)
      d_hocl = 6.4e-8 * tt / vis_h2so4
      k_hocl = 1.25e9 * ah * d_hocl * m_hcl
      c_hocl = 2009.0 * tt ** 0.5
      l_hocl = (d_hocl / k_hocl) ** 0.5
      f_hocl = 1 / TANH(radsa / l_hocl) - l_hocl / radsa
      if (ABS(f_hocl) < small_value .OR. ABS(f_hcl) < small_value) then
        gam3arr(jl,3) = 0.0
      else
        gam_hocl = (f_hocl * 4.0 * h_hocl * r_vol * tt / c_hocl) * (d_hocl * k_hocl) ** 0.5
        if (gam_hocl * f_hcl < small_value) then
          gam3arr(jl,3) = 0.0
        else
          gam_hocl_f_hcl = 1.0 + 1.0 / (gam_hocl * f_hcl)
          gam3arr(jl,3) = 1.0 / gam_hocl_f_hcl
          if (gam3arr(jl,3) < min_value) then
            gam3arr(jl,3) = min_value
          end if
        end if
      end if
      gam3arr(jl,8) = 1.0 / (1.0 / 0.8 + 1.0 / (EXP(29.2 - 0.4 * wt_h2so4) + 0.11))
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_shi_liquid_aerosol
  subroutine ukca_solidphase(n_points)
    use asad_mod, only : f, jpcspf, specf, sphno3
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: errcode
    integer :: js
    integer, save :: ihno3 = 0
    logical, save :: firstcall = .true.
    logical, save :: first_pass = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_SOLIDPHASE'
    integer :: idx

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first_pass) then
      if (firstcall) then
        !NO!$acc kernels
        !NO!$acc loop seq
        do js = 1, jpcspf, 1
          if (specf(js) == 'HONO2     ') then
            ihno3 = js
          end if
        enddo
        !NO!$acc end kernels
        if (ihno3 == 0) then
          errcode = 1
          cmessage(:) = 'Select HONO2 as advected tracer.'
          call ereport('SOLIDPHASE', errcode, cmessage)
        end if
        firstcall = .false.
      end if
      first_pass = .false.
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx = 1, n_points, 1
      f(idx,ihno3) = f(idx,ihno3) + sphno3(idx)
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_solidphase
  subroutine ukca_calckpsc(sasa, t, th2o, thcl, tcnit, tn2o5, thocl, akpsc1, akpsc2, akpsc3, akpsc4, akpsc5, lpsa, lphocl, lppsc, lpsimp, kchmlev, kstart, kend, dt)
    use asad_mod, only : fpsc1, sh2o, shno3
    use ukca_config_constants_mod, only : avogadro, boltzmann
    use ukca_constants, only : pi
    real, parameter :: u = 1.66056e-27
    real, parameter :: cpk1 = 1.0e+2
    real, parameter :: gam1a = 0.3
    real, parameter :: gam1b = 0.006
    real, parameter :: gam1c = 0.0006
    real, parameter :: gam1d = 0.003
    real, parameter :: gam1e = 0.3
    real, parameter :: gam2a = 0.3
    real, parameter :: gam2b = 0.3
    real, parameter :: gam2c = 0.03
    real, parameter :: gam2d = 0.03
    real, parameter :: gam2e = 0.3
    real, parameter :: gam3c = 0.1
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    logical, intent(in) :: lpsa
    logical, intent(in) :: lphocl
    logical, intent(in) :: lppsc
    logical, intent(in) :: lpsimp
    integer, intent(in) :: kchmlev
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    real, intent(in) :: dt
    real, dimension(kchmlev), intent(in) :: t
    real, dimension(kchmlev), intent(in) :: sasa
    real, dimension(kchmlev), intent(in) :: thcl
    real, dimension(kchmlev), intent(in) :: th2o
    real, dimension(kchmlev), intent(in) :: tcnit
    real, dimension(kchmlev), intent(in) :: tn2o5
    real, dimension(kchmlev), intent(in) :: thocl
    real, dimension(kchmlev), intent(out) :: akpsc1
    real, dimension(kchmlev), intent(out) :: akpsc2
    real, dimension(kchmlev), intent(out) :: akpsc3
    real, dimension(kchmlev), intent(out) :: akpsc4
    real, dimension(kchmlev), intent(out) :: akpsc5
    integer :: jl
    real :: zfcnit
    real :: zfn2o5
    real :: zfhocl
    real :: vol
    real :: order2
    real :: zrate
    real :: zdhcl
    real :: zfact
    real, dimension(kchmlev) :: gam3b
    real, dimension(kchmlev) :: hshcl
    real, dimension(kchmlev) :: hhocl
    real, dimension(kchmlev) :: ccnit
    real, dimension(kchmlev) :: cn2o5
    real, dimension(kchmlev) :: chocl
    real, dimension(kchmlev) :: psc1sa
    real, dimension(kchmlev) :: psc2sa
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CALCKPSC'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx = kstart, kend, 1
      akpsc1(idx) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = kstart, kend, 1
      akpsc2(idx_1) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_2 = kstart, kend, 1
      akpsc3(idx_2) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_3 = kstart, kend, 1
      akpsc4(idx_3) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_4 = kstart, kend, 1
      akpsc5(idx_4) = 0.0
    enddo
    !NO!$acc end kernels
    if (lppsc .OR. lpsa) then
      !NO!$acc kernels
      zfcnit = (8.0 * boltzmann) / (pi * 97.5 * u)
      zfn2o5 = (8.0 * boltzmann) / (pi * 108.0 * u)
      zfhocl = (8.0 * boltzmann) / (pi * 52.5 * u)
      !NO!$acc loop gang vector independent
      do idx_5 = kstart, kend, 1
        ccnit(idx_5) = 0.25 * SQRT(zfcnit * t(idx_5))
      enddo
      !NO!$acc loop gang vector independent
      do idx_6 = kstart, kend, 1
        cn2o5(idx_6) = 0.25 * SQRT(zfn2o5 * t(idx_6))
      enddo
      !NO!$acc loop gang vector independent
      do idx_7 = kstart, kend, 1
        chocl(idx_7) = 0.25 * SQRT(zfhocl * t(idx_7))
      enddo
      !NO!$acc end kernels
      if (lppsc) then
        if (lpsimp) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_8 = kstart, kend, 1
            akpsc1(idx_8) = 4.6e-5 * fpsc1(idx_8)
          enddo
          !NO!$acc loop gang vector independent
          do idx_9 = kstart, kend, 1
            akpsc2(idx_9) = 4.6e-5 * fpsc1(idx_9)
          enddo
          !NO!$acc loop gang vector independent
          do idx_10 = kstart, kend, 1
            akpsc3(idx_10) = 4.6e-5 * fpsc1(idx_10)
          enddo
          !NO!$acc end kernels
        else
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do idx_11 = kstart, kend, 1
            psc1sa(idx_11) = 1.85 * 63.0 * u * 3.0e5 * shno3(idx_11) / (rho1 * rad1)
          enddo
          !NO!$acc loop gang vector independent
          do idx_12 = kstart, kend, 1
            psc2sa(idx_12) = 18.0 * u * 3.0e5 * sh2o(idx_12) / (rho2 * rad2)
          enddo
          !NO!$acc loop gang vector independent
          do idx_13 = kstart, kend, 1
            akpsc1(idx_13) = ccnit(idx_13) * (psc1sa(idx_13) * gam1a + psc2sa(idx_13) * gam2a)
          enddo
          !NO!$acc loop gang vector independent
          do idx_14 = kstart, kend, 1
            akpsc2(idx_14) = ccnit(idx_14) * (psc1sa(idx_14) * gam1b + psc2sa(idx_14) * gam2b)
          enddo
          !NO!$acc loop gang vector independent
          do idx_15 = kstart, kend, 1
            akpsc3(idx_15) = cn2o5(idx_15) * (psc1sa(idx_15) * gam1c + psc2sa(idx_15) * gam2c)
          enddo
          !NO!$acc loop gang vector independent
          do idx_16 = kstart, kend, 1
            akpsc4(idx_16) = cn2o5(idx_16) * (psc1sa(idx_16) * gam1d + psc2sa(idx_16) * gam2d)
          enddo
          !NO!$acc loop gang vector independent
          do idx_17 = kstart, kend, 1
            akpsc5(idx_17) = chocl(idx_17) * (psc1sa(idx_17) * gam1e + psc2sa(idx_17) * gam2e)
          enddo
          !NO!$acc end kernels
        end if
      end if
      if (lpsa) then
        call ukca_eqcomp(t, th2o, kstart, kend, kchmlev, lphocl, gam3b, hshcl, hhocl)
        !NO!$acc kernels
        !NO!$acc loop gang vector independent
        do idx_18 = kstart, kend, 1
          akpsc2(idx_18) = akpsc2(idx_18) + ccnit(idx_18) * 100.0 * sasa(idx_18) * gam3b(idx_18)
        enddo
        !NO!$acc loop gang vector independent
        do idx_19 = kstart, kend, 1
          akpsc3(idx_19) = akpsc3(idx_19) + cn2o5(idx_19) * 100.0 * sasa(idx_19) * gam3c
        enddo
        !NO!$acc end kernels
        if (lphocl) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent
          do jl = kstart, kend, 1
            vol = sasa(jl) * radsa / 3.0
            order2 = 1.0e6 * cpk1 * vol * avogadro * hshcl(jl) * hhocl(jl) * (boltzmann * t(jl)) ** 2
            akpsc5(jl) = akpsc5(jl) + order2 * thcl(jl)
          enddo
          !NO!$acc end kernels
        end if
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent
      do jl = kstart, kend, 1
        zrate = MAX(1.0, dt * (akpsc1(jl) * tcnit(jl) + akpsc4(jl) * tn2o5(jl) + akpsc5(jl) * thocl(jl)))
        zdhcl = MIN(zrate, thcl(jl))
        zfact = zdhcl / zrate
        akpsc1(jl) = zfact * akpsc1(jl)
        akpsc4(jl) = zfact * akpsc4(jl)
        akpsc5(jl) = zfact * akpsc5(jl)
      enddo
      !NO!$acc end kernels
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_calckpsc
  subroutine ukca_eqcomp(t, th2o, kstart, kend, kchmlev, lphocl, rpcncl, hshcl, hhocl)
    use ukca_config_constants_mod, only : boltzmann
    integer, parameter :: noph2o = 16
    integer, parameter :: notemp = 28
    integer, parameter :: nocomp = 4
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kchmlev
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    logical, intent(in) :: lphocl
    real, dimension(kchmlev), intent(in) :: t
    real, dimension(kchmlev), intent(in) :: th2o
    real, dimension(kchmlev), intent(out) :: hshcl
    real, dimension(kchmlev), intent(out) :: rpcncl
    real, dimension(kchmlev), intent(out) :: hhocl
    integer :: jxp1
    integer :: jyp1
    integer :: jx
    integer :: jy
    integer :: ierx
    integer :: iery
    integer :: i
    integer :: jl
    real :: sxy
    real :: sx1y
    real :: sx1y1
    real :: sxy1
    real :: ta
    real :: tb
    real :: tt
    real :: ua
    real :: ub
    real :: u
    real :: a
    real :: b
    real :: ajx
    real :: ajxp1
    real :: bjx
    real :: bjxp1
    real :: ph2o
    real, dimension(kchmlev) :: comp
    real, dimension(nocomp), save :: compindx
    real, dimension(nocomp), save :: aindex
    real, dimension(nocomp), save :: bindex
    real, dimension(noph2o), save :: pindx
    real, dimension(notemp), save :: tindx
    real, dimension(noph2o,notemp), save :: c
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_EQCOMP'
    DATA compindx / 35.0, 40.0, 50.0, 60.0 /
    DATA aindex / -11.35, -11.95, -13.80, -14.68 /
    DATA bindex / 6.10E+03, 5.88E+03, 5.53E+03, 4.92E+03 /
    DATA pindx / 0.1000E-05, 0.1000E-04, 0.5000E-04, 0.1000E-03, 0.1500E-03, 0.2000E-03, 0.3000E-03, 0.5000E-03, 0.6000E-03, 0.8000E-03, 0.1000E-02, 0.1200E-02, 0.1500E-02, 0.2000E-02, 0.3000E-02, 0.1000E-01 /
    DATA tindx / 0.1750E+03, 0.1800E+03, 0.1850E+03, 0.1900E+03, 0.1950E+03, 0.2000E+03, 0.2050E+03, 0.2100E+03, 0.2150E+03, 0.2200E+03, 0.2250E+03, 0.2300E+03, 0.2350E+03, 0.2400E+03, 0.2450E+03, 0.2500E+03, 0.2550E+03, 0.2600E+03, 0.2650E+03, 0.2700E+03, 0.2750E+03, 0.2800E+03, 0.2850E+03, 0.2900E+03, 0.2950E+03, 0.3000E+03, 0.3050E+03, 0.3100E+03 /
    DATA (c(i, 1), i = 1, 16) / 0.9400E+02, 0.7800E+02, 0.5957E+02, 0.4363E+02, 0.3430E+02, 0.2768E+02, 0.1836E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01 /
    DATA (c(i, 2), i = 1, 16) / 0.9400E+02, 0.7900E+02, 0.6152E+02, 0.4723E+02, 0.3887E+02, 0.3294E+02, 0.2458E+02, 0.1272E+02, 0.4000E+02, 0.4000E+02, 0.1142E+02, 0.4000E+02, 0.1341E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01 /
    DATA (c(i, 3), i = 1, 16) / 0.9400E+02, 0.8000E+02, 0.6348E+02, 0.5084E+02, 0.4344E+02, 0.3820E+02, 0.3080E+02, 0.1929E+02, 0.1196E+02, 0.1681E+02, 0.1728E+02, 0.1300E+02, 0.1928E+02, 0.1442E+02, 0.4000E+02, 0.9000E+01 /
    DATA (c(i, 4), i = 1, 16) / 0.9400E+02, 0.8100E+02, 0.6543E+02, 0.5444E+02, 0.4801E+02, 0.4345E+02, 0.3702E+02, 0.2585E+02, 0.1538E+02, 0.2541E+02, 0.2315E+02, 0.1806E+02, 0.2515E+02, 0.1988E+02, 0.1615E+02, 0.9000E+01 /
    DATA (c(i, 5), i = 1, 16) / 0.9400E+02, 0.8200E+02, 0.6935E+02, 0.6165E+02, 0.5715E+02, 0.5396E+02, 0.4946E+02, 0.4226E+02, 0.3935E+02, 0.3402E+02, 0.2902E+02, 0.2313E+02, 0.3102E+02, 0.2535E+02, 0.2334E+02, 0.9000E+01 /
    DATA (c(i, 6), i = 1, 16) / 0.9400E+02, 0.8300E+02, 0.7125E+02, 0.6594E+02, 0.6283E+02, 0.6062E+02, 0.5751E+02, 0.5278E+02, 0.5073E+02, 0.4693E+02, 0.4369E+02, 0.4086E+02, 0.3689E+02, 0.3082E+02, 0.3052E+02, 0.9000E+01 /
    DATA (c(i, 7), i = 1, 16) / 0.9400E+02, 0.8367E+02, 0.7395E+02, 0.6976E+02, 0.6731E+02, 0.6557E+02, 0.6312E+02, 0.5955E+02, 0.5811E+02, 0.5561E+02, 0.5344E+02, 0.5144E+02, 0.4863E+02, 0.4449E+02, 0.3771E+02, 0.9000E+01 /
    DATA (c(i, 8), i = 1, 16) / 0.9400E+02, 0.8420E+02, 0.7626E+02, 0.7284E+02, 0.7084E+02, 0.6942E+02, 0.6742E+02, 0.6455E+02, 0.6341E+02, 0.6147E+02, 0.5983E+02, 0.5838E+02, 0.5646E+02, 0.5369E+02, 0.4849E+02, 0.1209E+02 /
    DATA (c(i, 9), i = 1, 16) / 0.9400E+02, 0.8519E+02, 0.7841E+02, 0.7548E+02, 0.7377E+02, 0.7256E+02, 0.7085E+02, 0.6845E+02, 0.6752E+02, 0.6594E+02, 0.6462E+02, 0.6347E+02, 0.6196E+02, 0.5983E+02, 0.5640E+02, 0.3239E+02 /
    DATA (c(i, 10), i = 1, 16) / 0.9400E+02, 0.8603E+02, 0.8020E+02, 0.7768E+02, 0.7621E+02, 0.7517E+02, 0.7370E+02, 0.7163E+02, 0.7083E+02, 0.6949E+02, 0.6839E+02, 0.6743E+02, 0.6619E+02, 0.6447E+02, 0.6175E+02, 0.4271E+02 /
    DATA (c(i, 11), i = 1, 16) / 0.9424E+02, 0.8691E+02, 0.8179E+02, 0.7959E+02, 0.7830E+02, 0.7738E+02, 0.7609E+02, 0.7429E+02, 0.7360E+02, 0.7244E+02, 0.7148E+02, 0.7066E+02, 0.6960E+02, 0.6815E+02, 0.6589E+02, 0.5007E+02 /
    DATA (c(i, 12), i = 1, 16) / 0.9433E+02, 0.8780E+02, 0.8323E+02, 0.8127E+02, 0.8012E+02, 0.7930E+02, 0.7815E+02, 0.7656E+02, 0.7595E+02, 0.7493E+02, 0.7410E+02, 0.7338E+02, 0.7245E+02, 0.7119E+02, 0.6925E+02, 0.5567E+02 /
    DATA (c(i, 13), i = 1, 16) / 0.9445E+02, 0.8860E+02, 0.8451E+02, 0.8275E+02, 0.8172E+02, 0.8099E+02, 0.7996E+02, 0.7853E+02, 0.7798E+02, 0.7708E+02, 0.7633E+02, 0.7570E+02, 0.7489E+02, 0.7377E+02, 0.7207E+02, 0.6017E+02 /
    DATA (c(i, 14), i = 1, 16) / 0.9478E+02, 0.8945E+02, 0.8571E+02, 0.8411E+02, 0.8317E+02, 0.8250E+02, 0.8156E+02, 0.8027E+02, 0.7977E+02, 0.7896E+02, 0.7829E+02, 0.7772E+02, 0.7699E+02, 0.7600E+02, 0.7449E+02, 0.6392E+02 /
    DATA (c(i, 15), i = 1, 16) / 0.9568E+02, 0.9057E+02, 0.8700E+02, 0.8546E+02, 0.8456E+02, 0.8392E+02, 0.8302E+02, 0.8183E+02, 0.8138E+02, 0.8063E+02, 0.8002E+02, 0.7951E+02, 0.7885E+02, 0.7795E+02, 0.7659E+02, 0.6707E+02 /
    DATA (c(i, 16), i = 1, 16) / 0.9695E+02, 0.9190E+02, 0.8836E+02, 0.8684E+02, 0.8595E+02, 0.8532E+02, 0.8443E+02, 0.8327E+02, 0.8284E+02, 0.8215E+02, 0.8158E+02, 0.8111E+02, 0.8050E+02, 0.7969E+02, 0.7845E+02, 0.6977E+02 /
    DATA (c(i, 17), i = 1, 16) / 0.9900E+02, 0.9374E+02, 0.9000E+02, 0.8840E+02, 0.8746E+02, 0.8679E+02, 0.8585E+02, 0.8467E+02, 0.8425E+02, 0.8357E+02, 0.8303E+02, 0.8258E+02, 0.8202E+02, 0.8126E+02, 0.8012E+02, 0.7214E+02 /
    DATA (c(i, 18), i = 1, 16) / 0.9900E+02, 0.9563E+02, 0.9170E+02, 0.9001E+02, 0.8902E+02, 0.8832E+02, 0.8733E+02, 0.8610E+02, 0.8566E+02, 0.8497E+02, 0.8444E+02, 0.8399E+02, 0.8344E+02, 0.8272E+02, 0.8164E+02, 0.7408E+02 /
    DATA (c(i, 19), i = 1, 16) / 0.9900E+02, 0.9753E+02, 0.9341E+02, 0.9163E+02, 0.9059E+02, 0.8985E+02, 0.8881E+02, 0.8753E+02, 0.8707E+02, 0.8637E+02, 0.8585E+02, 0.8540E+02, 0.8486E+02, 0.8418E+02, 0.8316E+02, 0.7602E+02 /
    DATA (c(i, 20), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9511E+02, 0.9324E+02, 0.9215E+02, 0.9138E+02, 0.9029E+02, 0.8896E+02, 0.8848E+02, 0.8777E+02, 0.8726E+02, 0.8681E+02, 0.8628E+02, 0.8564E+02, 0.8468E+02, 0.7796E+02 /
    DATA (c(i, 21), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9681E+02, 0.9486E+02, 0.9372E+02, 0.9291E+02, 0.9177E+02, 0.9039E+02, 0.8989E+02, 0.8917E+02, 0.8867E+02, 0.8822E+02, 0.8770E+02, 0.8710E+02, 0.8620E+02, 0.7990E+02 /
    DATA (c(i, 22), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9851E+02, 0.9647E+02, 0.9528E+02, 0.9444E+02, 0.9325E+02, 0.9182E+02, 0.9130E+02, 0.9057E+02, 0.9008E+02, 0.8963E+02, 0.8912E+02, 0.8856E+02, 0.8772E+02, 0.8184E+02 /
    DATA (c(i, 23), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9809E+02, 0.9685E+02, 0.9597E+02, 0.9473E+02, 0.9325E+02, 0.9271E+02, 0.9197E+02, 0.9149E+02, 0.9104E+02, 0.9054E+02, 0.9002E+02, 0.8924E+02, 0.8378E+02 /
    DATA (c(i, 24), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9842E+02, 0.9750E+02, 0.9621E+02, 0.9468E+02, 0.9412E+02, 0.9337E+02, 0.9290E+02, 0.9245E+02, 0.9196E+02, 0.9148E+02, 0.9076E+02, 0.8572E+02 /
    DATA (c(i, 25), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9769E+02, 0.9611E+02, 0.9553E+02, 0.9477E+02, 0.9431E+02, 0.9386E+02, 0.9338E+02, 0.9294E+02, 0.9228E+02, 0.8766E+02 /
    DATA (c(i, 26), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9754E+02, 0.9694E+02, 0.9617E+02, 0.9572E+02, 0.9527E+02, 0.9480E+02, 0.9440E+02, 0.9380E+02, 0.8960E+02 /
    DATA (c(i, 27), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9897E+02, 0.9835E+02, 0.9757E+02, 0.9713E+02, 0.9668E+02, 0.9622E+02, 0.9586E+02, 0.9532E+02, 0.9154E+02 /
    DATA (c(i, 28), i = 1, 16) / 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9897E+02, 0.9854E+02, 0.9809E+02, 0.9764E+02, 0.9732E+02, 0.9684E+02, 0.9348E+02 /
    integer :: idx
    integer :: idx_1
    integer :: idx_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx = LBOUND(hshcl, dim=1), UBOUND(hshcl, dim=1), 1
      hshcl(idx) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = LBOUND(rpcncl, dim=1), UBOUND(rpcncl, dim=1), 1
      rpcncl(idx_1) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_2 = LBOUND(hhocl, dim=1), UBOUND(hhocl, dim=1), 1
      hhocl(idx_2) = 0.0
    enddo
    !NO!$acc end kernels
    do jl = kstart, kend, 1
      ph2o = 1.0e4 * th2o(jl) * boltzmann * t(jl)
      call ukca_position(pindx, noph2o, ph2o, jx, ierx, 0)
      call ukca_position(tindx, notemp, t(jl), jy, iery, 0)
      if (jy == 0) then
        comp(jl) = 9.0
      else
        if (iery == 1) then
          comp(jl) = 95.0
        else
          if (jx == 0) then
            comp(jl) = 95.0
          else
            if (ierx == 1) then
              comp(jl) = 9.0
            else
              jxp1 = jx + 1
              jyp1 = jy + 1
              sxy = c(jx,jy)
              sx1y = c(jxp1,jy)
              sx1y1 = c(jxp1,jyp1)
              sxy1 = c(jx,jyp1)
              ta = ph2o - pindx(jx)
              tb = pindx(jxp1) - pindx(jx)
              tt = ta / tb
              ua = t(jl) - tindx(jy)
              ub = tindx(jyp1) - tindx(jy)
              u = ua / ub
              comp(jl) = (1.0 - tt) * (1.0 - u) * sxy + tt * (1.0 - u) * sx1y + tt * u * sx1y1 + (1.0 - tt) * u * sxy1
              if (comp(jl) < 9.0) then
                comp(jl) = 9.0
              end if
              if (comp(jl) > 95.0) then
                comp(jl) = 95.0
              end if
            end if
          end if
        end if
      end if
      rpcncl(jl) = 10.0 ** (1.87 - 0.074 * comp(jl))
      if (rpcncl(jl) > 1.0) then
        rpcncl(jl) = 1.0
      end if
    enddo
    if (lphocl) then
      do jl = kstart, kend, 1
        call ukca_position(compindx, nocomp, comp(jl), jx, ierx, 0)
        if (jx == 0) then
          a = -11.35
          b = 6.10e+03
        else
          if (ierx == 1) then
            a = -14.68
            b = 4.92e+03
          else
            jxp1 = jx + 1
            ajx = aindex(jx)
            ajxp1 = aindex(jxp1)
            bjx = bindex(jx)
            bjxp1 = bindex(jxp1)
            ta = comp(jl) - compindx(jx)
            tb = compindx(jxp1) - compindx(jx)
            tt = ta / tb
            a = (1.0 - tt) * ajx + tt * ajxp1
            b = (1.0 - tt) * bjx + tt * bjxp1
          end if
        end if
        hshcl(jl) = EXP(a + b / t(jl)) / 101.3
        hhocl(jl) = EXP(0.71 + 1633.0 / t(jl)) / 101.3
      enddo
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_eqcomp
  subroutine ukca_position(xc, n, x, jx, ier, iorder)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n
    integer, intent(in) :: iorder
    real, intent(in) :: x
    real, dimension(n), intent(in) :: xc
    integer, intent(out) :: ier
    integer, intent(out) :: jx
    integer :: i
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_POSITION'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    ier = 1
    if (iorder == 0) then
      if (x < xc(1)) then
        jx = 0
        ier = 0
      else
        do i = 1, n, 1
          if (x < xc(i)) then
            ier = 0
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
        jx = i - 1
      end if
    else
      if (x > xc(1)) then
        jx = 0
        ier = 0
      else
        do i = 1, n, 1
          if (x > xc(i)) then
            ier = 0
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Exit_Stmt
            EXIT
          end if
        enddo
        jx = i - 1
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_position
  subroutine ukca_pscpres(t, p, tnd, th2o, thno3, kstart, kend, kchmlev, have_nat, sph2o)
    use asad_mod, only : fpsc1, fpsc2, sh2o, shno3
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    integer, intent(in) :: kchmlev
    real, dimension(kchmlev), intent(in) :: t
    real, dimension(kchmlev), intent(in) :: p
    real, dimension(kchmlev), intent(in) :: tnd
    logical, dimension(kchmlev), intent(in) :: have_nat
    real, dimension(kchmlev), intent(inout) :: th2o
    real, dimension(kchmlev), intent(inout) :: thno3
    real, dimension(kchmlev), intent(inout) :: sph2o
    integer :: jl
    real :: zh2ot
    real :: ztpsc
    real :: zmt
    real :: zbt
    real :: zhno3eq
    real :: zh2oeq
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_PSCPRES'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent
    do idx = kstart, kend, 1
      shno3(idx) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_1 = kstart, kend, 1
      sh2o(idx_1) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_2 = kstart, kend, 1
      fpsc1(idx_2) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_3 = kstart, kend, 1
      fpsc2(idx_3) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do jl = kstart, kend, 1
      zh2ot = 0.75 * p(jl) * th2o(jl) / tnd(jl)
      zh2ot = MAX(zh2ot, 1.0e-12)
      ztpsc = t(jl)
      zmt = -2.7836 - 0.00088 * ztpsc
      zbt = 38.9855 - 11397.0 / ztpsc + 0.009179 * ztpsc
      zhno3eq = 10.0 ** (zmt * LOG10(zh2ot) + zbt)
      zhno3eq = zhno3eq * 133.3 * tnd(jl) / (100.0 * p(jl))
      if (thno3(jl) > zhno3eq .AND. have_nat(jl)) then
        fpsc1(jl) = 1.0
        shno3(jl) = thno3(jl) - zhno3eq
        thno3(jl) = zhno3eq
      else
        fpsc1(jl) = 0.0
      end if
      if (.true.) then
        if (sph2o(jl) > 0.0) then
          sh2o(jl) = sph2o(jl) * tnd(jl)
          fpsc2(jl) = 1.0
        else
          fpsc2(jl) = 0.0
          sh2o(jl) = 0.0
        end if
      else
        zh2oeq = 610.78 * EXP(21.875 * (t(jl) - 273.16) / (t(jl) - 7.66))
        zh2oeq = zh2oeq * tnd(jl) / (100.0 * p(jl))
        if (th2o(jl) > zh2oeq) then
          fpsc2(jl) = 1.0
          sh2o(jl) = th2o(jl) - zh2oeq
          th2o(jl) = zh2oeq
        else
          fpsc2(jl) = 0.0
        end if
      end if
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_pscpres

end module ukca_hetero_mod
