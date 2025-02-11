module get_molmass_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'GET_MOLMASS_MOD'
  public

  contains
  function get_molmass(species_name)
    use ukca_constants, only : m_benzene, m_br, m_brcl, m_bro, m_brono2, m_c, m_c2h2, m_c2h4, m_c2h6, m_c3h6, m_c3h8, m_c4h10, m_ccl4, m_cf2cl2, m_cf2clbr, m_cf2clcfcl2, m_cfcl3, m_ch2br2, m_ch4, m_chf2cl, m_cl, m_cl2o2, m_clo, m_clono2, m_co, m_cs2, m_dms, m_etcho, m_etoh, m_h2, m_h2s, m_hbr, m_hcho, m_hcl, m_hcooh, m_hobr, m_hoch2cho, m_hocl, m_isop, m_me2co, m_me2s, m_mebr, m_meccl3, m_mecho, m_mecl, m_meco2h, m_mek, m_meoh, m_monoterp, m_n2o, m_nh3, m_no, m_no2, m_oclo, m_ocs, m_oxylene, m_so2, m_so4, m_tbut2ene, m_toluene
    use ereport_mod, only : ereport
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    CHARACTER(LEN = 10), INTENT(IN) :: species_name
    real :: get_molmass
    integer :: ierr
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'GET_MOLMASS'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    ierr = 0
    if (species_name(:) == 'NO2       ' .OR. species_name(:) == 'NOx       ') then
      get_molmass = m_no2
    else
      if (species_name(:) == 'NO        ') then
        get_molmass = m_no
      else
        if (species_name(:) == 'H2        ' .OR. species_name(:) == 'no_H2     ') then
          get_molmass = m_h2
        else
          if (species_name(:) == 'NH3       ') then
            get_molmass = m_nh3
          else
            if (species_name(:) == 'N2O       ' .OR. species_name(:) == 'no_N2O    ') then
              get_molmass = m_n2o
            else
              if (species_name(:) == 'BC_fossil ' .OR. (species_name(:) == 'BC_biofuel' .OR. (species_name(:) == 'BC_biomass' .OR. (species_name(:) == 'OM_fossil ' .OR. (species_name(:) == 'OM_biofuel' .OR. species_name(:) == 'OM_biomass'))))) then
                get_molmass = m_c
              else
                if (species_name(:) == 'SO2       ' .OR. (species_name(:) == 'SO2_low   ' .OR. (species_name(:) == 'SO2_high  ' .OR. species_name(:) == 'SO2_nat   '))) then
                  get_molmass = m_so2
                else
                  if (species_name(:) == 'DMS       ') then
                    get_molmass = m_dms
                  else
                    if (species_name(:) == 'Me2S      ') then
                      get_molmass = m_me2s
                    else
                      if (species_name(:) == 'SO4       ') then
                        get_molmass = m_so4
                      else
                        if (species_name(:) == 'COS       ') then
                          get_molmass = m_ocs
                        else
                          if (species_name(:) == 'H2S       ') then
                            get_molmass = m_h2s
                          else
                            if (species_name(:) == 'CS2       ') then
                              get_molmass = m_cs2
                            else
                              if (species_name(:) == 'CO        ') then
                                get_molmass = m_co
                              else
                                if (species_name(:) == 'CH4       ') then
                                  get_molmass = m_ch4
                                else
                                  if (species_name(:) == 'C2H6      ') then
                                    get_molmass = m_c2h6
                                  else
                                    if (species_name(:) == 'C2H4      ') then
                                      get_molmass = m_c2h4
                                    else
                                      if (species_name(:) == 'C3H8      ') then
                                        get_molmass = m_c3h8
                                      else
                                        if (species_name(:) == 'C3H6      ') then
                                          get_molmass = m_c3h6
                                        else
                                          if (species_name(:) == 'C4H10     ') then
                                            get_molmass = m_c4h10
                                          else
                                            if (species_name(:) == 'MeOH      ' .OR. species_name(:) == 'CH3OH     ') then
                                              get_molmass = m_meoh
                                            else
                                              if (species_name(:) == 'Me2CO     ') then
                                                get_molmass = m_me2co
                                              else
                                                if (species_name(:) == 'HCHO      ') then
                                                  get_molmass = m_hcho
                                                else
                                                  if (species_name(:) == 'MeCHO     ') then
                                                    get_molmass = m_mecho
                                                  else
                                                    if (species_name(:) == 'C5H8      ') then
                                                      get_molmass = m_isop
                                                    else
                                                      if (species_name(:) == 'Monoterp  ' .OR. (species_name(:) == 'APINENE   ' .OR. species_name(:) == 'BPINENE   ')) then
                                                        get_molmass = m_monoterp
                                                      else
                                                        if (species_name(:) == 'TOLUENE   ') then
                                                          get_molmass = m_toluene
                                                        else
                                                          if (species_name(:) == 'oXYLENE   ') then
                                                            get_molmass = m_oxylene
                                                          else
                                                            if (species_name(:) == 'MeCl      ') then
                                                              get_molmass = m_mecl
                                                            else
                                                              if (species_name(:) == 'MeCCl3    ') then
                                                                get_molmass = m_meccl3
                                                              else
                                                                if (species_name(:) == 'CCl4      ') then
                                                                  get_molmass = m_ccl4
                                                                else
                                                                  if (species_name(:) == 'CF2Cl2    ') then
                                                                    get_molmass = m_cf2cl2
                                                                  else
                                                                    if (species_name(:) == 'CFCl3     ') then
                                                                      get_molmass = m_cfcl3
                                                                    else
                                                                      if (species_name(:) == 'CF2ClCFCl2') then
                                                                        get_molmass = m_cf2clcfcl2
                                                                      else
                                                                        if (species_name(:) == 'CF2ClBr   ') then
                                                                          get_molmass = m_cf2clbr
                                                                        else
                                                                          if (species_name(:) == 'CHF2Cl    ') then
                                                                            get_molmass = m_chf2cl
                                                                          else
                                                                            if (species_name(:) == 'Cl        ' .OR. species_name(:) == 'TOT_Cl    ') then
                                                                              get_molmass = m_cl
                                                                            else
                                                                              if (species_name(:) == 'ClO       ' .OR. species_name(:) == 'Clx       ') then
                                                                                get_molmass = m_clo
                                                                              else
                                                                                if (species_name(:) == 'Cl2O2     ') then
                                                                                  get_molmass = m_cl2o2
                                                                                else
                                                                                  if (species_name(:) == 'HCl       ') then
                                                                                    get_molmass = m_hcl
                                                                                  else
                                                                                    if (species_name(:) == 'HOCl      ') then
                                                                                      get_molmass = m_hocl
                                                                                    else
                                                                                      if (species_name(:) == 'ClONO2    ') then
                                                                                        get_molmass = m_clono2
                                                                                      else
                                                                                        if (species_name(:) == 'OClO      ') then
                                                                                          get_molmass = m_oclo
                                                                                        else
                                                                                          if (species_name(:) == 'MeBr      ') then
                                                                                            get_molmass = m_mebr
                                                                                          else
                                                                                            if (species_name(:) == 'CH2Br2    ') then
                                                                                              get_molmass = m_ch2br2
                                                                                            else
                                                                                              if (species_name(:) == 'Br        ' .OR. species_name(:) == 'TOT_Br    ') then
                                                                                                get_molmass = m_br
                                                                                              else
                                                                                                if (species_name(:) == 'BrO       ' .OR. species_name(:) == 'Brx       ') then
                                                                                                  get_molmass = m_bro
                                                                                                else
                                                                                                  if (species_name(:) == 'HBr       ') then
                                                                                                    get_molmass = m_hbr
                                                                                                  else
                                                                                                    if (species_name(:) == 'HOBr      ') then
                                                                                                      get_molmass = m_hobr
                                                                                                    else
                                                                                                      if (species_name(:) == 'BrONO2    ') then
                                                                                                        get_molmass = m_brono2
                                                                                                      else
                                                                                                        if (species_name(:) == 'BrCl      ') then
                                                                                                          get_molmass = m_brcl
                                                                                                        else
                                                                                                          if (species_name(:) == 'C2H2      ') then
                                                                                                            get_molmass = m_c2h2
                                                                                                          else
                                                                                                            if (species_name(:) == 'TBUT2ENE  ') then
                                                                                                              get_molmass = m_tbut2ene
                                                                                                            else
                                                                                                              if (species_name(:) == 'BENZENE   ') then
                                                                                                                get_molmass = m_benzene
                                                                                                              else
                                                                                                                if (species_name(:) == 'EtOH      ') then
                                                                                                                  get_molmass = m_etoh
                                                                                                                else
                                                                                                                  if (species_name(:) == 'EtCHO     ') then
                                                                                                                    get_molmass = m_etcho
                                                                                                                  else
                                                                                                                    if (species_name(:) == 'MEK       ') then
                                                                                                                      get_molmass = m_mek
                                                                                                                    else
                                                                                                                      if (species_name(:) == 'HCOOH     ') then
                                                                                                                        get_molmass = m_hcooh
                                                                                                                      else
                                                                                                                        if (species_name(:) == 'MeCO2H    ') then
                                                                                                                          get_molmass = m_meco2h
                                                                                                                        else
                                                                                                                          if (species_name(:) == 'HOCH2CHO  ') then
                                                                                                                            get_molmass = m_hoch2cho
                                                                                                                          else
                                                                                                                            if (species_name(:) == 'AGE       ') then
                                                                                                                              get_molmass = 1.0
                                                                                                                              ierr = -1
                                                                                                                            else
                                                                                                                              if (species_name(:) == 'AGE OF AIR') then
                                                                                                                                get_molmass = 1.0
                                                                                                                                ierr = -1
                                                                                                                              else
                                                                                                                                if (species_name(:) == 'PASSIVE O3') then
                                                                                                                                  get_molmass = 1.0
                                                                                                                                  ierr = -1
                                                                                                                                else
                                                                                                                                  if (species_name(:) == 'XXX       ') then
                                                                                                                                    get_molmass = 1.0
                                                                                                                                    ierr = -1
                                                                                                                                  else
                                                                                                                                    ierr = 1
                                                                                                                                  end if
                                                                                                                                end if
                                                                                                                              end if
                                                                                                                            end if
                                                                                                                          end if
                                                                                                                        end if
                                                                                                                      end if
                                                                                                                    end if
                                                                                                                  end if
                                                                                                                end if
                                                                                                              end if
                                                                                                            end if
                                                                                                          end if
                                                                                                        end if
                                                                                                      end if
                                                                                                    end if
                                                                                                  end if
                                                                                                end if
                                                                                              end if
                                                                                            end if
                                                                                          end if
                                                                                        end if
                                                                                      end if
                                                                                    end if
                                                                                  end if
                                                                                end if
                                                                              end if
                                                                            end if
                                                                          end if
                                                                        end if
                                                                      end if
                                                                    end if
                                                                  end if
                                                                end if
                                                              end if
                                                            end if
                                                          end if
                                                        end if
                                                      end if
                                                    end if
                                                  end if
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    if (ierr /= 0) then
      if (ierr == (-1)) then
        cmessage(:) = 'Molecular mass of ' // TRIM(species_name) // ' set to 1'
      else
        if (ierr == 1) then
          cmessage(:) = 'Add molecular mass of ' // TRIM(species_name) // ' to this routine and check if it is present in ukca_constants'
        end if
      end if
      call ereport('GET_MOLMASS', ierr, cmessage)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end function get_molmass

end module get_molmass_mod
