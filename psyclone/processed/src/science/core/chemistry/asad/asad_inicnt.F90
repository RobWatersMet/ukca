module asad_inicnt_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INICNT_MOD'
  public

  contains
  subroutine asad_inicnt(species, y_out, klen, nlev)
    use asad_mod, only : co2, f, jpro2, nlfro2, tnd, wp
    use ukca_config_specification_mod, only : ukca_config
    use ukca_constants, only : c_ho2, c_no3, c_o3, c_oh
    use ukca_environment_fields_mod, only : ho2_offline, no3_offline, o3_offline, oh_offline
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: klen
    integer, intent(in) :: nlev
    CHARACTER(LEN = 10), INTENT(IN) :: species
    real, dimension(klen), intent(out) :: y_out
    integer :: errcode
    integer :: row_length
    integer :: rows
    integer :: iro2
    integer :: j
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    real, dimension(klen) :: fro2
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INICNT'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    row_length = ukca_config%row_length
    rows = ukca_config%rows
    if (species(1:4) == 'CO2 ' .AND. ukca_config%l_chem_environ_co2_fld) then
      do idx = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
        y_out(idx) = co2(idx + (LBOUND(co2, dim=1) - LBOUND(y_out, dim=1)))
      enddo
    else
      if (species(1:4) == 'H2O ') then
        do idx_1 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
          y_out(idx_1) = wp(idx_1 + (LBOUND(wp, dim=1) - LBOUND(y_out, dim=1)))
        enddo
      else
        if (species(1:4) == 'RO2 ') then
          if (ukca_config%l_ukca_ro2_perm) then
            do idx_2 = LBOUND(fro2, dim=1), UBOUND(fro2, dim=1), 1
              fro2(idx_2) = 0.0
            enddo
            do j = 1, jpro2, 1
              iro2 = nlfro2(j)
              do idx_3 = LBOUND(fro2, dim=1), UBOUND(fro2, dim=1), 1
                fro2(idx_3) = fro2(idx_3) + f(idx_3 + (LBOUND(f, dim=1) - LBOUND(fro2, dim=1)),iro2)
              enddo
            enddo
            do idx_4 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
              y_out(idx_4) = fro2(idx_4) / tnd(idx_4 + (LBOUND(tnd, dim=1) - LBOUND(y_out, dim=1)))
            enddo
          else
            errcode = 126
            cmessage(:) = 'RO2 should only be a species if l_ukca_ro2_perm == T'
            call ereport('ASAD_INICNT', errcode, cmessage)
          end if
        else
          if (ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be) then
            if (species(1:4) == 'OH  ') then
              y_out(:) = RESHAPE(oh_offline(1:row_length,1:rows,nlev), [klen])
              y_out(:) = y_out(:) / c_oh
            else
              if (species(1:4) == 'O3  ') then
                y_out(:) = RESHAPE(o3_offline(1:row_length,1:rows,nlev), [klen])
                y_out(:) = y_out(:) / c_o3
              else
                if (species(1:4) == 'NO3 ') then
                  y_out(:) = RESHAPE(no3_offline(1:row_length,1:rows,nlev), [klen])
                  y_out(:) = y_out(:) / c_no3
                else
                  if (species(1:4) == 'HO2 ') then
                    y_out(:) = RESHAPE(ho2_offline(1:row_length,1:rows,nlev), [klen])
                    y_out(:) = y_out(:) / c_ho2
                  else
                    errcode = 125
                    cmessage(:) = ' Species ' // species // ' is not treated by this routine'
                    call ereport('ASAD_INICNT', errcode, cmessage)
                  end if
                end if
              end if
            end if
          else
            errcode = 124
            cmessage(:) = ' Species ' // species // ' not treated by this routine'
            call ereport('ASAD_INICNT', errcode, cmessage)
          end if
        end if
      end if
    end if
    do idx_5 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
      y_out(idx_5) = y_out(idx_5) * tnd(idx_5 + (LBOUND(tnd, dim=1) - LBOUND(y_out, dim=1)))
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_inicnt

end module asad_inicnt_mod
