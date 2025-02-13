module asad_fyinit_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_FYINIT_MOD'
  public

  contains
  subroutine asad_fyinit(ofirst, n_points, ix, jy, nlev)
    use asad_mod, only : ctype, f, ilcf, ilct, ilftr, ilss, jpcf, jpco, jpif, jpna, jpoo, jpsp, jpspec, jsro2, madvtr, moffam, nlmajmin, peps, speci, tnd, y
    use ukca_diurnal_oxidant, only : ukca_set_diurnal_ox, ukca_set_diurnal_ox_col
    use ukca_config_specification_mod, only : ukca_config
    use ukca_environment_fields_mod, only : atmospheric_ch4, atmospheric_co2, atmospheric_h2, atmospheric_n2, atmospheric_o2
    use ukca_constants, only : c_ch4, c_co2, c_h2, c_n2, c_o2
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    use asad_inicnt_mod, only : asad_inicnt
    use asad_inicnt_col_mod, only : asad_inicnt_col
    real, parameter :: fco2_default = 350.0e-6
    real, parameter :: fh2_default = 5.0e-7
    real, parameter :: fn2_default = 0.78084
    real, parameter :: fo2_default = 0.20945
    real, parameter :: fch4_default = 1.76e-6
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    logical, intent(in) :: ofirst
    integer, save :: iss
    integer, save :: ict
    integer, save :: iftr
    integer, save :: icf
    integer :: j
    integer :: jl
    integer :: js
    integer :: istart
    integer :: iend
    integer :: ifam
    integer :: itr
    integer :: errcode
    real :: peps10
    real :: fch4
    real :: fco2
    real :: fh2
    real :: fn2
    real :: fo2
    real :: zcnst
    logical, save :: gonce = .true.
    logical, save :: first_pass = .true.
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_FYINIT'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first_pass) then
      if (gonce) then
        !$acc kernels
        gonce = .false.
        iss = 0
        ict = 0
        icf = 0
        iftr = 0
        !$acc loop seq
        do js = 1, jpspec, 1
          if (ctype(js) == jpna(:)) then
            iss = iss + 1
            ilss(iss) = js
          end if
          if (ctype(js) == jpco(:)) then
            ict = ict + 1
            ilct(ict) = js
          end if
          if (ctype(js) == jpcf(:)) then
            icf = icf + 1
            ilcf(icf) = js
            if (speci(js) == 'RO2       ') then
              jsro2 = js
            end if
          end if
          if (ctype(js) == jpif(:) .OR. ctype(js) == jpsp(:) .OR. ctype(js) == jpoo(:)) then
            iftr = iftr + 1
            ilftr(iftr) = js
          end if
        enddo
        !$acc end kernels
      end if
      first_pass = .false.
    end if
    if (ofirst) then
      !$acc kernels
      peps10 = 10.0 * peps
      istart = nlmajmin(1)
      iend = nlmajmin(2)
      !$acc loop seq
      do j = istart, iend, 1
        js = nlmajmin(j)
        ifam = moffam(js)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          y(jl,js) = f(jl,ifam)
        enddo
      enddo
      istart = nlmajmin(3)
      iend = nlmajmin(4)
      !$acc loop seq
      do j = istart, iend, 1
        js = nlmajmin(j)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          y(jl,js) = peps10
        enddo
      enddo
      !$acc loop seq
      do j = 1, iss, 1
        js = ilss(j)
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          y(jl,js) = peps10
        enddo
      enddo
      fco2 = fco2_default
      fch4 = fch4_default
      fh2 = fh2_default
      fo2 = fo2_default
      fn2 = fn2_default
      !$acc end kernels
      if (ukca_config%l_chem_environ_co2_scalar) then
        fco2 = atmospheric_co2 / c_co2
        if (ukca_config%l_chem_environ_co2_fld) then
          fco2 = 0.0
        end if
      end if
      if (ukca_config%l_chem_environ_h2_scalar) then
        fh2 = atmospheric_h2 / c_h2
      end if
      if (ukca_config%l_chem_environ_n2_scalar) then
        fn2 = atmospheric_n2 / c_n2
      end if
      if (ukca_config%l_chem_environ_o2_scalar) then
        fo2 = atmospheric_o2 / c_o2
      end if
      if (ukca_config%l_chem_environ_ch4_scalar) then
        fch4 = atmospheric_ch4 / c_ch4
      end if
      do j = 1, ict, 1
        js = ilct(j)
        if (speci(js) == 'N2        ') then
          zcnst = fn2
        else
          if (speci(js) == 'O2        ') then
            zcnst = fo2
          else
            if (speci(js) == 'CO2       ') then
              zcnst = fco2
            else
              if (speci(js) == 'H2        ') then
                zcnst = fh2
              else
                if (speci(js) == 'CH4       ') then
                  zcnst = fch4
                else
                  cmessage(:) = ' Value not supplied for ' // speci(js) // ',setting value to zero'
                  errcode = -1
                  call ereport('ASAD_FYINIT', errcode, cmessage)
                  zcnst = 0.0
                end if
              end if
            end if
          end if
        end if
        !$acc kernels
        !$acc loop gang vector independent
        do jl = 1, n_points, 1
          y(jl,js) = zcnst * tnd(jl)
        enddo
        !$acc end kernels
      enddo
    end if
    !$acc kernels
    !$acc loop seq
    do j = 1, iftr, 1
      js = ilftr(j)
      itr = madvtr(js)
      !$acc loop gang vector independent
      do jl = 1, n_points, 1
        y(jl,js) = f(jl,itr)
      enddo
    enddo
    !$acc end kernels
    if (ofirst) then
      do j = 1, icf, 1
        js = ilcf(j)
        if (ukca_config%l_ukca_asad_columns) then
          call asad_inicnt_col(speci(js), y(1,js), n_points, ix, jy)
        else
          call asad_inicnt(speci(js), y(1,js), n_points, nlev)
        end if
      enddo
    end if
    if (ofirst .AND. (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be)) then
      do j = 1, icf, 1
        js = ilcf(j)
        if (ukca_config%l_ukca_asad_columns) then
          call ukca_set_diurnal_ox_col(speci(js), y(:,js), n_points, ix, jy)
        else
          call ukca_set_diurnal_ox(speci(js), y(:,js), n_points, nlev)
        end if
      enddo
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_fyinit

end module asad_fyinit_mod
