module asad_inicnt_col_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_INICNT_COL_MOD'
  public

  contains
  subroutine asad_inicnt_col(species, y_out, klen, ix, jy)
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
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    CHARACTER(LEN = 10), INTENT(IN) :: species
    real, dimension(klen), intent(out) :: y_out
    integer :: errcode
    integer :: iro2
    integer :: j
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    real, dimension(klen) :: fro2
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_INICNT_COL'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
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
              do idx_5 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
                y_out(idx_5) = oh_offline(ix,jy,idx_5 + (LBOUND(oh_offline, dim=3) - LBOUND(y_out, dim=1)))
              enddo
              y_out(:) = y_out(:) / c_oh
            else
              if (species(1:4) == 'O3  ') then
                do idx_6 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
                  y_out(idx_6) = o3_offline(ix,jy,idx_6 + (LBOUND(o3_offline, dim=3) - LBOUND(y_out, dim=1)))
                enddo
                y_out(:) = y_out(:) / c_o3
              else
                if (species(1:4) == 'NO3 ') then
                  do idx_7 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
                    y_out(idx_7) = no3_offline(ix,jy,idx_7 + (LBOUND(no3_offline, dim=3) - LBOUND(y_out, dim=1)))
                  enddo
                  y_out(:) = y_out(:) / c_no3
                else
                  if (species(1:4) == 'HO2 ') then
                    do idx_8 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
                      y_out(idx_8) = ho2_offline(ix,jy,idx_8 + (LBOUND(ho2_offline, dim=3) - LBOUND(y_out, dim=1)))
                    enddo
                    y_out(:) = y_out(:) / c_ho2
                  else
                    errcode = 125
                    cmessage(:) = ' Species ' // species // ' is not treated by this routine'
                    call ereport('ASAD_INICNT_COL', errcode, cmessage)
                  end if
                end if
              end if
            end if
          else
            errcode = 124
            cmessage(:) = ' Species ' // species // ' not treated by this routine'
            call ereport('ASAD_INICNT_COL', errcode, cmessage)
          end if
        end if
      end if
    end if
    do idx_9 = LBOUND(y_out, dim=1), UBOUND(y_out, dim=1), 1
      y_out(idx_9) = y_out(idx_9) * tnd(idx_9 + (LBOUND(tnd, dim=1) - LBOUND(y_out, dim=1)))
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine asad_inicnt_col

end module asad_inicnt_col_mod
