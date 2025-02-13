module ukca_aerod_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_AEROD_MOD'
  public

  contains
  subroutine ukca_aerod(row_length, rows, ntype, npft, nwater, t0, p0, hflx, u_s, canht, gsf, zbl, z0tile, resa, rb, so4_vd)
    use asad_mod, only : jpdd, ndepd, nldepd, speci
    use ukca_config_specification_mod, only : ukca_config
    use ukca_constants, only : m_air, m_alkaooh, m_anhy, m_arnoh14, m_arnoh17, m_aroh14, m_aroh17, m_aromooh, m_asoa, m_asvoc1, m_asvoc2, m_bsoa, m_bsvoc1, m_bsvoc2, m_buooh, m_carb10, m_carb11a, m_carb12, m_carb13, m_carb14, m_carb15, m_carb16, m_carb17, m_carb9, m_ccarb12, m_cri, m_dhcarb9, m_dhpcarb9, m_dhpr12ooh, m_dmso, m_etcho, m_etco3h, m_etoh, m_gly, m_h2o, m_h2so4, m_hbr, m_hcho, m_hcl, m_hmml, m_hobr, m_hoc2h4no3, m_hoc2h4ooh, m_hoch2cho, m_hoch2co3h, m_hocl, m_hono, m_hpucarb12, m_hucarb9, m_iepox, m_isooh, m_isosoa, m_isosvoc1, m_isosvoc2, m_macr, m_macrooh, m_mecho, m_mecoch2ooh, m_mekooh, m_meoh, m_mgly, m_monoterp, m_mpan, m_msa, m_mvkooh, m_nald, m_noa, m_nrn12ooh, m_nrn6ooh, m_nrn9ooh, m_nrtn28ooh, m_nrtx28ooh, m_nru12ooh, m_nru14ooh, m_nucarb12, m_onitu, m_pan, m_phan, m_ppan, m_proh, m_ra13no3, m_ra13ooh, m_ra16no3, m_ra16ooh, m_ra19no3, m_ra19ooh, m_rcooh25, m_rn12no3, m_rn12ooh, m_rn14ooh, m_rn15no3, m_rn15ooh, m_rn16ooh, m_rn17ooh, m_rn18no3, m_rn18ooh, m_rn19ooh, m_rn9no3, m_rn9ooh, m_rtn10ooh, m_rtn14ooh, m_rtn23no3, m_rtn23ooh, m_rtn24ooh, m_rtn25no3, m_rtn25ooh, m_rtn26ooh, m_rtn26pan, m_rtn28no3, m_rtn28ooh, m_rtx22no3, m_rtx22ooh, m_rtx24no3, m_rtx24ooh, m_rtx28no3, m_rtx28ooh, m_ru10no3, m_ru12no3, m_ru12pan, m_sec_org, m_sec_org_i, m_tncarb10, m_tncarb11, m_tncarb12, m_tncarb15, m_tncarb26, m_txcarb22, m_txcarb24, m_ucarb12, m_udcarb11, m_udcarb14, m_udcarb8
    use ukca_um_legacy_mod, only : gg=>g, vkman
    use ukca_config_constants_mod, only : rmol
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    use nvtx
    real, parameter :: r_null = 1.0e50
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: ntype
    integer, intent(in) :: npft
    integer, intent(in) :: nwater
    real, dimension(row_length,rows), intent(in) :: t0
    real, dimension(row_length,rows), intent(in) :: p0
    real, dimension(row_length,rows), intent(in) :: hflx
    real, dimension(row_length,rows), intent(in) :: u_s
    real, dimension(row_length,rows,npft), intent(in) :: canht
    real, dimension(row_length,rows,ntype), intent(in) :: gsf
    real, dimension(row_length,rows), intent(in) :: zbl
    real, dimension(row_length,rows,ntype), intent(inout) :: z0tile
    real, dimension(row_length,rows,ntype), intent(out) :: resa
    real, dimension(row_length,rows,jpdd), intent(out) :: rb
    real, dimension(row_length,rows), intent(out) :: so4_vd
    logical, save :: first = .true.
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    real, save :: cp = 1004.67
    real, save :: pr = 0.72
    real, save :: vair296 = 1830.0e-8
    real, save :: twoth = 2.0 / 3.0
    real, save :: zref = 50.0
    real, save :: d_h2o = 2.08e-5
    real :: b0
    real :: b1
    real :: sc_pr
    real, save :: reus
    real, allocatable, dimension(:), save :: d0
    real, dimension(row_length,rows) :: l
    real, dimension(row_length,rows) :: rho_air
    real, dimension(row_length,rows) :: kva
    real, dimension(row_length,rows) :: ustar
    real, dimension(row_length,rows,npft) :: d
    real, dimension(row_length,rows,ntype) :: z
    real, dimension(row_length,rows,ntype) :: psi
    real :: bl_l
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_AEROD'
    CHARACTER(LEN = errormessagelength) :: cmessage
    integer :: icode
    integer :: idx
    integer :: idx_1
    integer :: idx_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (first) then
      ALLOCATE(d0(1:jpdd))
      !$acc kernels
      !$acc loop gang vector independent
      do idx = LBOUND(d0, dim=1), UBOUND(d0, dim=1), 1
        d0(idx) = -1.0
      enddo
      !$acc end kernels
      !$acc kernels
      !$acc loop gang vector independent
      do j = 1, ndepd, 1
        if (speci(nldepd(j)) == 'O3        ' .OR. (speci(nldepd(j)) == 'NO2       ' .OR. (speci(nldepd(j)) == 'O3S       ' .OR. speci(nldepd(j)) == 'NO3       '))) then
          d0(j) = 1.4e-5
        else
          if (speci(nldepd(j)) == 'NO        ') then
            if (ukca_config%l_fix_improve_drydep) then
              d0(j) = 2.3e-5
            else
              d0(j) = 8.4e-5
            end if
          else
            if (speci(nldepd(j)) == 'HNO3      ' .OR. (speci(nldepd(j)) == 'HONO2     ' .OR. (speci(nldepd(j)) == 'ISON      ' .OR. (speci(nldepd(j)) == 'B2ndry    ' .OR. (speci(nldepd(j)) == 'A2ndry    ' .OR. (speci(nldepd(j)) == 'N2O5      ' .OR. (speci(nldepd(j)) == 'HO2NO2    ' .OR. (speci(nldepd(j)) == 'HNO4      ' .OR. speci(nldepd(j)) == 'RU14NO3   ')))))))) then
              d0(j) = 1.2e-5
            else
              if (speci(nldepd(j)) == 'H2O2      ' .OR. speci(nldepd(j)) == 'HOOH      ') then
                d0(j) = 1.46e-5
              else
                if (speci(nldepd(j)) == 'CH3OOH    ' .OR. (speci(nldepd(j)) == 'MeOOH     ' .OR. speci(nldepd(j)) == 'HCOOH     ')) then
                  d0(j) = 1.27e-5
                else
                  if (speci(nldepd(j)) == 'C2H5OOH   ' .OR. (speci(nldepd(j)) == 'EtOOH     ' .OR. (speci(nldepd(j)) == 'MeCO3H    ' .OR. (speci(nldepd(j)) == 'MeCO2H    ' .OR. (speci(nldepd(j)) == 'HACET     ' .OR. speci(nldepd(j)) == 'CARB7     '))))) then
                    d0(j) = 1.12e-5
                  else
                    if (speci(nldepd(j)) == 'n_C3H7OOH ' .OR. (speci(nldepd(j)) == 'i_C3H7OOH ' .OR. (speci(nldepd(j)) == 'n-PrOOH   ' .OR. (speci(nldepd(j)) == 'i-PrOOH   ' .OR. (speci(nldepd(j)) == 'PropeOOH  ' .OR. speci(nldepd(j)) == 'RN10OOH   '))))) then
                      d0(j) = 1.01e-5
                    else
                      if (speci(nldepd(j)) == 'MeCOCH2OOH') then
                        d0(j) = d_h2o * SQRT(m_h2o / m_mecoch2ooh)
                      else
                        if (speci(nldepd(j)) == 'ISOOH     ' .OR. speci(nldepd(j)) == 'RU14OOH   ') then
                          d0(j) = d_h2o * SQRT(m_h2o / m_isooh)
                        else
                          if (speci(nldepd(j)) == 'HONO      ') then
                            d0(j) = d_h2o * SQRT(m_h2o / m_hono)
                          else
                            if (speci(nldepd(j)) == 'MACROOH   ' .OR. speci(nldepd(j)) == 'RU10OOH   ') then
                              d0(j) = d_h2o * SQRT(m_h2o / m_macrooh)
                            else
                              if (speci(nldepd(j)) == 'MEKOOH    ' .OR. speci(nldepd(j)) == 'RN11OOH   ') then
                                d0(j) = d_h2o * SQRT(m_h2o / m_mekooh)
                              else
                                if (speci(nldepd(j)) == 'ALKAOOH   ' .OR. speci(nldepd(j)) == 'RN8OOH    ') then
                                  d0(j) = d_h2o * SQRT(m_h2o / m_alkaooh)
                                else
                                  if (speci(nldepd(j)) == 'AROMOOH   ') then
                                    d0(j) = d_h2o * SQRT(m_h2o / m_aromooh)
                                  else
                                    if (speci(nldepd(j)) == 'BSVOC1    ') then
                                      d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc1)
                                    else
                                      if (speci(nldepd(j)) == 'BSVOC2    ') then
                                        d0(j) = d_h2o * SQRT(m_h2o / m_bsvoc2)
                                      else
                                        if (speci(nldepd(j)) == 'ASVOC1    ') then
                                          d0(j) = d_h2o * SQRT(m_h2o / m_asvoc1)
                                        else
                                          if (speci(nldepd(j)) == 'ASVOC2    ') then
                                            d0(j) = d_h2o * SQRT(m_h2o / m_asvoc2)
                                          else
                                            if (speci(nldepd(j)) == 'ISOSVOC1  ') then
                                              d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc1)
                                            else
                                              if (speci(nldepd(j)) == 'ISOSVOC2  ') then
                                                d0(j) = d_h2o * SQRT(m_h2o / m_isosvoc2)
                                              else
                                                if (speci(nldepd(j)) == 'PAN       ') then
                                                  if (ukca_config%l_fix_improve_drydep) then
                                                    d0(j) = d_h2o * SQRT(m_h2o / m_pan)
                                                  else
                                                    d0(j) = 0.31e-5
                                                  end if
                                                else
                                                  if (speci(nldepd(j)) == 'PPAN      ') then
                                                    d0(j) = d_h2o * SQRT(m_h2o / m_ppan)
                                                  else
                                                    if (speci(nldepd(j)) == 'MPAN      ') then
                                                      d0(j) = d_h2o * SQRT(m_h2o / m_mpan)
                                                    else
                                                      if (speci(nldepd(j)) == 'ONITU     ') then
                                                        d0(j) = d_h2o * SQRT(m_h2o / m_onitu)
                                                      else
                                                        if (speci(nldepd(j)) == 'CO        ') then
                                                          d0(j) = 1.86e-5
                                                        else
                                                          if (speci(nldepd(j)) == 'CH4       ') then
                                                            if (ukca_config%l_fix_improve_drydep) then
                                                              d0(j) = 2.2e-5
                                                            else
                                                              d0(j) = 5.74e-5
                                                            end if
                                                          else
                                                            if (speci(nldepd(j)) == 'NH3       ') then
                                                              d0(j) = 2.08e-5
                                                            else
                                                              if (speci(nldepd(j)) == 'H2        ') then
                                                                d0(j) = 6.7e-5
                                                              else
                                                                if (speci(nldepd(j)) == 'SO2       ') then
                                                                  d0(j) = 1.2e-5
                                                                else
                                                                  if (speci(nldepd(j)) == 'DMSO      ') then
                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_dmso)
                                                                  else
                                                                    if (speci(nldepd(j)) == 'Sec_Org   ') then
                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_sec_org)
                                                                    else
                                                                      if (speci(nldepd(j)) == 'SEC_ORG_I ') then
                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_sec_org_i)
                                                                      else
                                                                        if (speci(nldepd(j)) == 'H2SO4     ') then
                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_h2so4)
                                                                        else
                                                                          if (speci(nldepd(j)) == 'MSA       ') then
                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_msa)
                                                                          else
                                                                            if (speci(nldepd(j)) == 'MVKOOH    ') then
                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_mvkooh)
                                                                            else
                                                                              if (speci(nldepd(j)) == 's-BuOOH   ' .OR. speci(nldepd(j)) == 'RN13OOH   ') then
                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_buooh)
                                                                              else
                                                                                if (speci(nldepd(j)) == 'ORGNIT    ') then
                                                                                  d0(j) = 1.0e-5
                                                                                else
                                                                                  if (speci(nldepd(j)) == 'BSOA      ') then
                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_bsoa)
                                                                                  else
                                                                                    if (speci(nldepd(j)) == 'ASOA      ') then
                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_asoa)
                                                                                    else
                                                                                      if (speci(nldepd(j)) == 'ISOSOA    ') then
                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_isosoa)
                                                                                      else
                                                                                        if (speci(nldepd(j)) == 'HCHO      ') then
                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_hcho)
                                                                                        else
                                                                                          if (speci(nldepd(j)) == 'MeCHO     ') then
                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_mecho)
                                                                                          else
                                                                                            if (speci(nldepd(j)) == 'EtCHO     ') then
                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_etcho)
                                                                                            else
                                                                                              if (speci(nldepd(j)) == 'MACR      ' .OR. speci(nldepd(j)) == 'UCARB10   ') then
                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_macr)
                                                                                              else
                                                                                                if (speci(nldepd(j)) == 'NALD      ') then
                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_nald)
                                                                                                else
                                                                                                  if (speci(nldepd(j)) == 'GLY       ' .OR. speci(nldepd(j)) == 'CARB3     ') then
                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_gly)
                                                                                                  else
                                                                                                    if (speci(nldepd(j)) == 'MGLY      ' .OR. speci(nldepd(j)) == 'CARB6     ') then
                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_mgly)
                                                                                                    else
                                                                                                      if (speci(nldepd(j)) == 'MeOH      ') then
                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_meoh)
                                                                                                      else
                                                                                                        if (speci(nldepd(j)) == 'Monoterp  ' .OR. (speci(nldepd(j)) == 'APINENE   ' .OR. speci(nldepd(j)) == 'BPINENE   ')) then
                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_monoterp)
                                                                                                        else
                                                                                                          if (speci(nldepd(j)) == 'HBr       ') then
                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_hbr)
                                                                                                          else
                                                                                                            if (speci(nldepd(j)) == 'HOBr      ') then
                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_hobr)
                                                                                                            else
                                                                                                              if (speci(nldepd(j)) == 'HCl       ') then
                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_hcl)
                                                                                                              else
                                                                                                                if (speci(nldepd(j)) == 'HOCl      ') then
                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_hocl)
                                                                                                                else
                                                                                                                  if (speci(nldepd(j)) == 'EtOH      ') then
                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_etoh)
                                                                                                                  else
                                                                                                                    if (speci(nldepd(j)) == 'i-PrOH    ' .OR. speci(nldepd(j)) == 'n-PrOH    ') then
                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_proh)
                                                                                                                    else
                                                                                                                      if (speci(nldepd(j)) == 'HOC2H4OOH ') then
                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_hoc2h4ooh)
                                                                                                                      else
                                                                                                                        if (speci(nldepd(j)) == 'HOCH2CHO  ') then
                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_hoch2cho)
                                                                                                                        else
                                                                                                                          if (speci(nldepd(j)) == 'EtCO3H    ') then
                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_etco3h)
                                                                                                                          else
                                                                                                                            if (speci(nldepd(j)) == 'HOCH2CO3H ') then
                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_hoch2co3h)
                                                                                                                            else
                                                                                                                              if (speci(nldepd(j)) == 'PHAN      ') then
                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_phan)
                                                                                                                              else
                                                                                                                                if (speci(nldepd(j)) == 'NOA       ') then
                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_noa)
                                                                                                                                else
                                                                                                                                  if (speci(nldepd(j)) == 'HOC2H4NO3 ') then
                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_hoc2h4no3)
                                                                                                                                  else
                                                                                                                                    if (speci(nldepd(j)) == 'RTX24NO3  ') then
                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_rtx24no3)
                                                                                                                                    else
                                                                                                                                      if (speci(nldepd(j)) == 'RN9NO3    ') then
                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_rn9no3)
                                                                                                                                      else
                                                                                                                                        if (speci(nldepd(j)) == 'RN12NO3   ') then
                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rn12no3)
                                                                                                                                        else
                                                                                                                                          if (speci(nldepd(j)) == 'RN15NO3   ') then
                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_rn15no3)
                                                                                                                                          else
                                                                                                                                            if (speci(nldepd(j)) == 'RN18NO3   ') then
                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_rn18no3)
                                                                                                                                            else
                                                                                                                                              if (speci(nldepd(j)) == 'RTN28NO3  ') then
                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_rtn28no3)
                                                                                                                                              else
                                                                                                                                                if (speci(nldepd(j)) == 'RTN25NO3  ') then
                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_rtn25no3)
                                                                                                                                                else
                                                                                                                                                  if (speci(nldepd(j)) == 'RTN23NO3  ') then
                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_rtn23no3)
                                                                                                                                                  else
                                                                                                                                                    if (speci(nldepd(j)) == 'RTX28NO3  ') then
                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_rtx28no3)
                                                                                                                                                    else
                                                                                                                                                      if (speci(nldepd(j)) == 'RTX22NO3  ') then
                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_rtx22no3)
                                                                                                                                                      else
                                                                                                                                                        if (speci(nldepd(j)) == 'RN16OOH   ') then
                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rn16ooh)
                                                                                                                                                        else
                                                                                                                                                          if (speci(nldepd(j)) == 'RN19OOH   ') then
                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_rn19ooh)
                                                                                                                                                          else
                                                                                                                                                            if (speci(nldepd(j)) == 'RN14OOH   ') then
                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_rn14ooh)
                                                                                                                                                            else
                                                                                                                                                              if (speci(nldepd(j)) == 'RN17OOH   ') then
                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_rn17ooh)
                                                                                                                                                              else
                                                                                                                                                                if (speci(nldepd(j)) == 'NRU14OOH  ') then
                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_nru14ooh)
                                                                                                                                                                else
                                                                                                                                                                  if (speci(nldepd(j)) == 'NRU12OOH  ') then
                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_nru12ooh)
                                                                                                                                                                  else
                                                                                                                                                                    if (speci(nldepd(j)) == 'RN9OOH    ') then
                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_rn9ooh)
                                                                                                                                                                    else
                                                                                                                                                                      if (speci(nldepd(j)) == 'RN12OOH   ') then
                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_rn12ooh)
                                                                                                                                                                      else
                                                                                                                                                                        if (speci(nldepd(j)) == 'RN15OOH   ') then
                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rn15ooh)
                                                                                                                                                                        else
                                                                                                                                                                          if (speci(nldepd(j)) == 'RN18OOH   ') then
                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_rn18ooh)
                                                                                                                                                                          else
                                                                                                                                                                            if (speci(nldepd(j)) == 'NRN6OOH   ') then
                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_nrn6ooh)
                                                                                                                                                                            else
                                                                                                                                                                              if (speci(nldepd(j)) == 'NRN9OOH   ') then
                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_nrn9ooh)
                                                                                                                                                                              else
                                                                                                                                                                                if (speci(nldepd(j)) == 'NRN12OOH  ') then
                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_nrn12ooh)
                                                                                                                                                                                else
                                                                                                                                                                                  if (speci(nldepd(j)) == 'RA13OOH   ') then
                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_ra13ooh)
                                                                                                                                                                                  else
                                                                                                                                                                                    if (speci(nldepd(j)) == 'RA16OOH   ') then
                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_ra16ooh)
                                                                                                                                                                                    else
                                                                                                                                                                                      if (speci(nldepd(j)) == 'RA19OOH   ') then
                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_ra19ooh)
                                                                                                                                                                                      else
                                                                                                                                                                                        if (speci(nldepd(j)) == 'RTN28OOH  ') then
                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rtn28ooh)
                                                                                                                                                                                        else
                                                                                                                                                                                          if (speci(nldepd(j)) == 'NRTN28OOH ') then
                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_nrtn28ooh)
                                                                                                                                                                                          else
                                                                                                                                                                                            if (speci(nldepd(j)) == 'RTN26OOH  ') then
                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_rtn26ooh)
                                                                                                                                                                                            else
                                                                                                                                                                                              if (speci(nldepd(j)) == 'RTN25OOH  ') then
                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_rtn25ooh)
                                                                                                                                                                                              else
                                                                                                                                                                                                if (speci(nldepd(j)) == 'RTN24OOH  ') then
                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_rtn24ooh)
                                                                                                                                                                                                else
                                                                                                                                                                                                  if (speci(nldepd(j)) == 'RTN23OOH  ') then
                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_rtn23ooh)
                                                                                                                                                                                                  else
                                                                                                                                                                                                    if (speci(nldepd(j)) == 'RTN14OOH  ') then
                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_rtn14ooh)
                                                                                                                                                                                                    else
                                                                                                                                                                                                      if (speci(nldepd(j)) == 'RTN10OOH  ') then
                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_rtn10ooh)
                                                                                                                                                                                                      else
                                                                                                                                                                                                        if (speci(nldepd(j)) == 'RTX28OOH  ') then
                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rtx28ooh)
                                                                                                                                                                                                        else
                                                                                                                                                                                                          if (speci(nldepd(j)) == 'RTX24OOH  ') then
                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_rtx24ooh)
                                                                                                                                                                                                          else
                                                                                                                                                                                                            if (speci(nldepd(j)) == 'RTX22OOH  ') then
                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_rtx22ooh)
                                                                                                                                                                                                            else
                                                                                                                                                                                                              if (speci(nldepd(j)) == 'NRTX28OOH ') then
                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_nrtx28ooh)
                                                                                                                                                                                                              else
                                                                                                                                                                                                                if (speci(nldepd(j)) == 'CARB14    ') then
                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_carb14)
                                                                                                                                                                                                                else
                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'CARB17    ') then
                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_carb17)
                                                                                                                                                                                                                  else
                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'CARB11A   ') then
                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_carb11a)
                                                                                                                                                                                                                    else
                                                                                                                                                                                                                      if (speci(nldepd(j)) == 'CARB10    ') then
                                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_carb10)
                                                                                                                                                                                                                      else
                                                                                                                                                                                                                        if (speci(nldepd(j)) == 'CARB13    ') then
                                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_carb13)
                                                                                                                                                                                                                        else
                                                                                                                                                                                                                          if (speci(nldepd(j)) == 'CARB16    ') then
                                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_carb16)
                                                                                                                                                                                                                          else
                                                                                                                                                                                                                            if (speci(nldepd(j)) == 'CARB9     ') then
                                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_carb9)
                                                                                                                                                                                                                            else
                                                                                                                                                                                                                              if (speci(nldepd(j)) == 'CARB12    ') then
                                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_carb12)
                                                                                                                                                                                                                              else
                                                                                                                                                                                                                                if (speci(nldepd(j)) == 'CARB15    ') then
                                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_carb15)
                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'UCARB12   ') then
                                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_ucarb12)
                                                                                                                                                                                                                                  else
                                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'NUCARB12  ') then
                                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_nucarb12)
                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                      if (speci(nldepd(j)) == 'UDCARB8   ') then
                                                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_udcarb8)
                                                                                                                                                                                                                                      else
                                                                                                                                                                                                                                        if (speci(nldepd(j)) == 'UDCARB11  ') then
                                                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_udcarb11)
                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                          if (speci(nldepd(j)) == 'UDCARB14  ') then
                                                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_udcarb14)
                                                                                                                                                                                                                                          else
                                                                                                                                                                                                                                            if (speci(nldepd(j)) == 'TNCARB26  ') then
                                                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_tncarb26)
                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                              if (speci(nldepd(j)) == 'TNCARB10  ') then
                                                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_tncarb10)
                                                                                                                                                                                                                                              else
                                                                                                                                                                                                                                                if (speci(nldepd(j)) == 'TNCARB12  ') then
                                                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_tncarb12)
                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'TNCARB11  ') then
                                                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_tncarb11)
                                                                                                                                                                                                                                                  else
                                                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'CCARB12   ') then
                                                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_ccarb12)
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                      if (speci(nldepd(j)) == 'TNCARB15  ') then
                                                                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_tncarb15)
                                                                                                                                                                                                                                                      else
                                                                                                                                                                                                                                                        if (speci(nldepd(j)) == 'RCOOH25   ') then
                                                                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_rcooh25)
                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                          if (speci(nldepd(j)) == 'TXCARB24  ') then
                                                                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_txcarb24)
                                                                                                                                                                                                                                                          else
                                                                                                                                                                                                                                                            if (speci(nldepd(j)) == 'TXCARB22  ') then
                                                                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_txcarb22)
                                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                              if (speci(nldepd(j)) == 'RU12PAN   ') then
                                                                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_ru12pan)
                                                                                                                                                                                                                                                              else
                                                                                                                                                                                                                                                                if (speci(nldepd(j)) == 'RTN26PAN  ') then
                                                                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_rtn26pan)
                                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'AROH14    ') then
                                                                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_aroh14)
                                                                                                                                                                                                                                                                  else
                                                                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'AROH17    ') then
                                                                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_aroh17)
                                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                                      if (speci(nldepd(j)) == 'ARNOH14   ') then
                                                                                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_arnoh14)
                                                                                                                                                                                                                                                                      else
                                                                                                                                                                                                                                                                        if (speci(nldepd(j)) == 'ARNOH17   ') then
                                                                                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_arnoh17)
                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                          if (speci(nldepd(j)) == 'ANHY      ') then
                                                                                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_anhy)
                                                                                                                                                                                                                                                                          else
                                                                                                                                                                                                                                                                            if (speci(nldepd(j)) == 'IEPOX     ') then
                                                                                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_iepox)
                                                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                                              if (speci(nldepd(j)) == 'HMML      ') then
                                                                                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_hmml)
                                                                                                                                                                                                                                                                              else
                                                                                                                                                                                                                                                                                if (speci(nldepd(j)) == 'HUCARB9   ') then
                                                                                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_hucarb9)
                                                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'HPUCARB12 ') then
                                                                                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_hpucarb12)
                                                                                                                                                                                                                                                                                  else
                                                                                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'DHPCARB9  ') then
                                                                                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_dhpcarb9)
                                                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                                                      if (speci(nldepd(j)) == 'DHPR12OOH ') then
                                                                                                                                                                                                                                                                                        d0(j) = d_h2o * SQRT(m_h2o / m_dhpr12ooh)
                                                                                                                                                                                                                                                                                      else
                                                                                                                                                                                                                                                                                        if (speci(nldepd(j)) == 'DHCARB9   ') then
                                                                                                                                                                                                                                                                                          d0(j) = d_h2o * SQRT(m_h2o / m_dhcarb9)
                                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                                          if (speci(nldepd(j)) == 'RU12NO3   ') then
                                                                                                                                                                                                                                                                                            d0(j) = d_h2o * SQRT(m_h2o / m_ru12no3)
                                                                                                                                                                                                                                                                                          else
                                                                                                                                                                                                                                                                                            if (speci(nldepd(j)) == 'RU10NO3  ') then
                                                                                                                                                                                                                                                                                              d0(j) = d_h2o * SQRT(m_h2o / m_ru10no3)
                                                                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                                                              if (speci(nldepd(j)) == 'RU12OOH   ') then
                                                                                                                                                                                                                                                                                                d0(j) = d_h2o * SQRT(m_h2o / m_cri)
                                                                                                                                                                                                                                                                                              else
                                                                                                                                                                                                                                                                                                if (speci(nldepd(j)) == 'RA13NO3   ') then
                                                                                                                                                                                                                                                                                                  d0(j) = d_h2o * SQRT(m_h2o / m_ra13no3)
                                                                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                                                                  if (speci(nldepd(j)) == 'RA16NO3   ') then
                                                                                                                                                                                                                                                                                                    d0(j) = d_h2o * SQRT(m_h2o / m_ra16no3)
                                                                                                                                                                                                                                                                                                  else
                                                                                                                                                                                                                                                                                                    if (speci(nldepd(j)) == 'RA19NO3   ') then
                                                                                                                                                                                                                                                                                                      d0(j) = d_h2o * SQRT(m_h2o / m_ra19no3)
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
      enddo
      !$acc end kernels
      do j = 1, ndepd, 1
        if (d0(j) < 0.0) then
          cmessage(:) = 'Warning: No dry deposition for ' // speci(nldepd(j)) // ' will be calculated.'
          icode = -1
          call ereport(routinename, icode, cmessage)
        end if
      enddo
      reus = -cp / (vkman * gg)
      first = .false.
    end if
    !$acc kernels
    !$acc loop gang vector independent collapse(2)
    do k = 1, rows, 1
      do i = 1, row_length, 1
        rho_air(i,k) = m_air * p0(i,k) / (rmol * t0(i,k))
        ustar(i,k) = u_s(i,k)
        if (ustar(i,k) <= 0.0) then
          ustar(i,k) = 1.0e-2
        end if
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do n = 1, npft, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          if (canht(i,k,n) > 0.0) then
            d(i,k,n) = canht(i,k,n) * 0.70
          else
            d(i,k,n) = 0.0
          end if
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do n = 1, ntype, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          if (z0tile(i,k,n) <= 0.0) then
            z0tile(i,k,n) = 0.001
          end if
          z(i,k,n) = zref
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do n = 1, npft, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          z(i,k,n) = zref - d(i,k,n)
        enddo
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do k = 1, rows, 1
      do i = 1, row_length, 1
        kva(i,k) = (vair296 - 4.83e-8 * (296.0 - t0(i,k))) / rho_air(i,k)
        if (gsf(i,k,nwater) > 0.0) then
          z0tile(i,k,nwater) = kva(i,k) / (9.1 * ustar(i,k)) + 0.016 * ustar(i,k) * ustar(i,k) / gg
        end if
      enddo
    enddo
    !$acc loop gang vector independent collapse(2)
    do k = 1, rows, 1
      do i = 1, row_length, 1
        if (hflx(i,k) /= 0.0) then
          l(i,k) = reus * t0(i,k) * rho_air(i,k) * ustar(i,k) ** 3 / hflx(i,k)
        else
          l(i,k) = 10000.0
        end if
      enddo
    enddo
    !$acc loop seq
    do k = 1, rows, 1
      !$acc loop seq
      do i = 1, row_length, 1
        if (l(i,k) > 0.0) then
          !$acc loop gang vector independent
          do n = 1, ntype, 1
            psi(i,k,n) = -5.0 * (z(i,k,n) - z0tile(i,k,n)) / l(i,k)
          enddo
        else
          !$acc loop gang vector independent
          do n = 1, ntype, 1
            b1 = (1.0 - 16.0 * (z(i,k,n) / l(i,k))) ** 0.5
            b0 = (1.0 - 16.0 * (z0tile(i,k,n) / l(i,k))) ** 0.5
            psi(i,k,n) = 2.0 * LOG((1.0 + b1) / (1.0 + b0))
          enddo
        end if
      enddo
    enddo
    !$acc loop gang vector independent collapse(3)
    do n = 1, ntype, 1
      do k = 1, rows, 1
        do i = 1, row_length, 1
          resa(i,k,n) = (LOG(z(i,k,n) / z0tile(i,k,n)) - psi(i,k,n)) / (vkman * ustar(i,k))
        enddo
      enddo
    enddo
    !$acc loop seq
    do j = 1, ndepd, 1
      if (speci(nldepd(j)) == 'ORGNIT    ' .OR. speci(nldepd(j)) == 'ASOA      ' .OR. speci(nldepd(j)) == 'BSOA      ' .OR. speci(nldepd(j)) == 'ISOSOA    ') then
        !$acc loop gang vector independent collapse(2)
        do k = 1, rows, 1
          do i = 1, row_length, 1
            rb(i,k,j) = 1.0
          enddo
        enddo
      else
        if (d0(j) > 0.0) then
          !$acc loop gang vector independent collapse(2)
          do k = 1, rows, 1
            do i = 1, row_length, 1
              sc_pr = kva(i,k) / (pr * d0(j))
              rb(i,k,j) = sc_pr ** twoth / (vkman * ustar(i,k))
            enddo
          enddo
        else
          !$acc loop gang vector independent collapse(2)
          do k = 1, rows, 1
            do i = 1, row_length, 1
              rb(i,k,j) = 0.0
            enddo
          enddo
        end if
      end if
    enddo
    !$acc loop gang vector independent collapse(2)
    do idx_1 = 1, rows, 1
      do idx_2 = 1, row_length, 1
        so4_vd(idx_2,idx_1) = 1.0 / r_null
      enddo
    enddo
    !$acc end kernels
    do j = 1, ndepd, 1
      if (speci(nldepd(j)) == 'ORGNIT    ' .OR. speci(nldepd(j)) == 'ASOA      ' .OR. speci(nldepd(j)) == 'BSOA      ' .OR. speci(nldepd(j)) == 'ISOSOA    ') then
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do k = 1, rows, 1
          do i = 1, row_length, 1
            so4_vd(i,k) = 0.0
            bl_l = zbl(i,k) / l(i,k)
            if (l(i,k) >= 0.0) then
              so4_vd(i,k) = 0.002 * ustar(i,k)
            end if
            if (bl_l < (-30.0)) then
              so4_vd(i,k) = 0.0009 * ustar(i,k) * (-bl_l) ** twoth
            end if
            if (bl_l >= (-30.0) .AND. l(i,k) < 0.0) then
              so4_vd(i,k) = 0.002 * ustar(i,k) * (1.0 + (-300.0 / l(i,k)) ** twoth)
            end if
          enddo
        enddo
        !$acc end kernels
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Exit_Stmt
        EXIT
      end if
    enddo
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_aerod

end module ukca_aerod_mod
