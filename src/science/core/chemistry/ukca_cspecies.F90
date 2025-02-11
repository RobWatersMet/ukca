module ukca_cspecies
  use yomhook, only : dr_hook, lhook
  use parkind1, only : jpim, jprb
  implicit none
  integer, save, public :: n_ox
  integer, save, public :: n_o3
  integer, save, public :: n_o3s
  integer, save, public :: n_nox
  integer, save, public :: n_no
  integer, save, public :: n_no2
  integer, save, public :: n_no3
  integer, save, public :: n_n2o5
  integer, save, public :: n_ho2no2
  integer, save, public :: n_hono2
  integer, save, public :: n_pan
  integer, save, public :: n_ison
  integer, save, public :: n_orgnit
  integer, save, public :: n_rnc2h4
  integer, save, public :: n_rnc3h6
  integer, save, public :: n_h2o2
  integer, save, public :: n_ch4
  integer, save, public :: n_sx
  integer, save, public :: n_h2
  integer, save, public :: n_h2o
  integer, save, public :: n_cl
  integer, save, public :: n_ox1
  integer, save, public :: n_sx1
  integer, save, public :: n_dms
  integer, save, public :: n_so2
  integer, save, public :: n_so3
  integer, save, public :: n_h2so4
  integer, save, public :: n_sec_org
  integer, save, public :: n_sec_org_i
  integer, save, public :: n_cfcl3
  integer, save, public :: n_cf2cl2
  integer, save, public :: n_bro
  integer, save, public :: n_hcl
  integer, save, public :: n_o1d
  integer, save, public :: n_passive
  integer, save, public :: n_hcho
  integer, save, public :: n_c2h6
  integer, save, public :: n_mecho
  integer, save, public :: n_c3h8
  integer, save, public :: n_me2co
  integer, save, public :: n_c5h8
  integer, save, public :: n_mgly
  integer, save, public :: n_c4h10
  integer, save, public :: n_mek
  integer, save, public :: n_c2h4
  integer, save, public :: n_c3h6
  integer, save, public :: n_oxylene
  integer, save, public :: n_ch3oh
  integer, save, public :: n_gly
  integer, save, public :: n_mvk
  integer, save, public :: n_toluene
  integer, save, public :: n_co
  integer, save, public :: n_n2o
  integer, save, public :: n_mebr
  integer, save, public :: n_nh3
  integer, save, public :: nn_o3
  integer, save, public :: nn_o3s
  integer, save, public :: nn_oh
  integer, save, public :: nn_ho2
  integer, save, public :: nn_h2o2
  integer, save, public :: nn_no
  integer, save, public :: nn_no2
  integer, save, public :: nn_o1d
  integer, save, public :: nn_o3p
  integer, save, public :: nn_meoo
  integer, save, public :: nn_meco3
  integer, save, public :: nn_etoo
  integer, save, public :: nn_etco3
  integer, save, public :: nn_nproo
  integer, save, public :: nn_nprooh
  integer, save, public :: nn_iproo
  integer, save, public :: nn_iprooh
  integer, save, public :: nn_mecoch2oo
  integer, save, public :: nn_no3
  integer, save, public :: nn_n2o5
  integer, save, public :: nn_ho2no2
  integer, save, public :: nn_hono2
  integer, save, public :: nn_ch4
  integer, save, public :: nn_so2
  integer, save, public :: nn_so3
  integer, save, public :: nn_h2so4
  integer, save, public :: nn_ohs
  integer, save, public :: nn_ho2s
  integer, save, public :: nn_o1ds
  integer, save, public :: nn_o3ps
  integer, save, public :: nn_n
  integer, save, public :: nn_cl
  integer, save, public :: nn_clo
  integer, save, public :: nn_hcl
  integer, save, public :: nn_cl2o2
  integer, save, public :: nn_meo
  integer, save, public :: nn_bro
  integer, save, public :: nn_br
  integer, save, public :: nn_buoo
  integer, save, public :: nn_meko2
  integer, save, public :: nn_hoc2h4o2
  integer, save, public :: nn_hoc3h6o2
  integer, save, public :: nn_oxyl1
  integer, save, public :: nn_memald1
  integer, save, public :: nn_hoipo2
  integer, save, public :: nn_homvko2
  integer, save, public :: nn_tolp1
  integer, save, public :: nn_cfcl3
  integer, save, public :: nn_cf2cl2
  integer, save, public :: nn_nh3
  real, allocatable, dimension(:), save, public :: c_species
  real, allocatable, dimension(:), save, public :: c_na_species
  character(len = *), parameter, private :: ModuleName = 'UKCA_CSPECIES'
  private

  public :: ukca_calc_cspecies

  contains
  subroutine ukca_calc_cspecies()
    use ukca_constants, only : c_alkaooh, c_anhy, c_arnoh14, c_arnoh17, c_aroh14, c_aroh17, c_benzene, c_br, c_brcl, c_bro, c_brono2, c_buoo, c_buooh, c_c2h2, c_c2h4, c_c2h6, c_c3h6, c_c3h8, c_c4h10, c_c5h8, c_carb10, c_carb11a, c_carb12, c_carb13, c_carb14, c_carb15, c_carb16, c_carb17, c_carb9, c_ccarb12, c_ccl4, c_cf2cl2, c_cf2clbr, c_cf2clcfcl2, c_cf3br, c_cfcl3, c_ch2br2, c_ch3oh, c_ch3s, c_ch3sch2oo, c_ch3so, c_ch3so2, c_ch3so3, c_ch4, c_chf2cl, c_cl, c_cl2o2, c_clo, c_clono2, c_co, c_cos, c_cri, c_cs2, c_dhcarb9, c_dhpcarb9, c_dhpr12o2, c_dhpr12ooh, c_dms, c_dmso, c_etcho, c_etco3, c_etco3h, c_etoh, c_etono2, c_etoo, c_etooh, c_gly, c_h, c_h2, c_h2o, c_h2o2, c_h2s, c_h2so4, c_hacet, c_hbr, c_hcho, c_hcl, c_hcooh, c_hmml, c_ho2, c_ho2no2, c_hobr, c_hoc2h4no3, c_hoc2h4o2, c_hoc2h4ooh, c_hoc3h6o2, c_hoch2ch2o2, c_hoch2cho, c_hoch2co3, c_hoch2co3h, c_hocl, c_hoipo2, c_homvko2, c_hono, c_hono2, c_hpucarb12, c_hucarb9, c_iepox, c_iso2, c_ison, c_isooh, c_maco3, c_macr, c_macro2, c_macrooh, c_me2co, c_me2s, c_mebr, c_meccl3, c_mecho, c_mecl, c_meco2h, c_meco3, c_meco3h, c_mecoch2oo, c_mecoch2ooh, c_mek, c_meko2, c_mekooh, c_memald, c_memald1, c_meo2no2, c_meoh, c_meono2, c_meoo, c_meooh, c_mgly, c_monoterp, c_mpan, c_msa, c_msia, c_mvk, c_mvkooh, c_n, c_n2o, c_n2o5, c_nald, c_nh3, c_no, c_no2, c_no3, c_noa, c_nrn12ooh, c_nrn6ooh, c_nrn9ooh, c_nrtn28ooh, c_nrtx28ooh, c_nru12ooh, c_nru14ooh, c_nucarb12, c_o1d, c_o3, c_o3p, c_oclo, c_oh, c_orgnit, c_oxyl1, c_oxylene, c_pan, c_phan, c_ppan, c_proh, c_prono2, c_proo, c_prooh, c_ra13no3, c_ra13ooh, c_ra16no3, c_ra16ooh, c_ra19no3, c_ra19ooh, c_rcooh25, c_rn12no3, c_rn12ooh, c_rn13no3, c_rn14ooh, c_rn15no3, c_rn15ooh, c_rn16no3, c_rn16ooh, c_rn17ooh, c_rn18no3, c_rn18ooh, c_rn19no3, c_rn19ooh, c_rn9no3, c_rn9ooh, c_rnc2h4, c_rnc3h6, c_rtn10ooh, c_rtn14ooh, c_rtn23no3, c_rtn23ooh, c_rtn24ooh, c_rtn25no3, c_rtn25ooh, c_rtn26ooh, c_rtn26pan, c_rtn28no3, c_rtn28ooh, c_rtx22no3, c_rtx22ooh, c_rtx24no3, c_rtx24ooh, c_rtx28no3, c_rtx28ooh, c_ru10ao2, c_ru10no3, c_ru12no3, c_ru12pan, c_ru14no3, c_sec_org, c_sec_org_i, c_so2, c_so3, c_tbut2ene, c_tncarb10, c_tncarb11, c_tncarb12, c_tncarb15, c_tncarb26, c_tolp1, c_toluene, c_txcarb22, c_txcarb24, c_ucarb12, c_udcarb11, c_udcarb14, c_udcarb8
    use asad_mod, only : advt, jpctr, jpspec, nadvt, speci
    use ukca_config_specification_mod, only : ukca_config
    use ereport_mod, only : ereport
    use umPrintMgr, only : printstatus, prstatus_oper, ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer :: errcode
    integer :: i
    integer :: m
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CALC_CSPECIES'
    integer :: idx
    integer :: idx_1

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    !$acc loop gang vector independent
    do idx = LBOUND(c_species, dim=1), UBOUND(c_species, dim=1), 1
      c_species(idx) = 0.0
    enddo
    !$acc end kernels
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (advt == 'Ox        ') c_species = c_o3
    WHERE (advt == 'O3        ') c_species = c_o3
    WHERE (advt == 'NO        ') c_species = c_no
    WHERE (advt == 'NO2       ') c_species = c_no2
    WHERE (advt == 'NO3       ') c_species = c_no3
    WHERE (advt == 'NOx       ') c_species = c_no2
    WHERE (advt == 'N2O5      ') c_species = c_n2o5
    WHERE (advt == 'HO2NO2    ') c_species = c_ho2no2
    WHERE (advt == 'HONO2     ') c_species = c_hono2
    WHERE (advt == 'H2O2      ') c_species = c_h2o2
    WHERE (advt == 'CH4       ') c_species = c_ch4
    WHERE (advt == 'CO        ') c_species = c_co
    WHERE (advt == 'HCHO      ') c_species = c_hcho
    WHERE (advt == 'MeOOH     ') c_species = c_meooh
    WHERE (advt == 'HONO      ') c_species = c_hono
    WHERE (advt == 'C2H6      ') c_species = c_c2h6
    WHERE (advt == 'EtOOH     ') c_species = c_etooh
    WHERE (advt == 'MeCHO     ') c_species = c_mecho
    WHERE (advt == 'PAN       ') c_species = c_pan
    WHERE (advt == 'C3H8      ') c_species = c_c3h8
    WHERE (advt == 'i-PrOOH   ') c_species = c_prooh
    WHERE (advt == 'n-PrOOH   ') c_species = c_prooh
    WHERE (advt == 'EtCHO     ') c_species = c_etcho
    WHERE (advt == 'Me2CO     ') c_species = c_me2co
    WHERE (advt == 'MeCOCH2OOH') c_species = c_mecoch2ooh
    WHERE (advt == 'PPAN      ') c_species = c_ppan
    WHERE (advt == 'MeONO2    ') c_species = c_meono2
    WHERE (advt == 'Sx        ') c_species = c_o3
    WHERE (advt == 'O3S       ') c_species = c_o3
    WHERE (advt == 'HOx       ') c_species = c_ho2
    WHERE (advt == 'N2O       ') c_species = c_n2o
    WHERE (advt == 'CFCl3     ') c_species = c_cfcl3
    WHERE (advt == 'CF2Cl2    ') c_species = c_cf2cl2
    WHERE (advt == 'H2O       ') c_species = c_h2o
    WHERE (advt == 'ClONO2    ') c_species = c_clono2
    WHERE (advt == 'Clx       ') c_species = c_clo
    WHERE (advt == 'ClO       ') c_species = c_clo
    WHERE (advt == 'Cl2O2     ') c_species = c_cl2o2
    WHERE (advt == 'Cl        ') c_species = c_cl
    WHERE (advt == 'HCl       ') c_species = c_hcl
    WHERE (advt == 'HOCl      ') c_species = c_hocl
    WHERE (advt == 'OClO      ') c_species = c_oclo
    WHERE (advt == 'Br        ') c_species = c_br
    WHERE (advt == 'BrO       ') c_species = c_bro
    WHERE (advt == 'Brx       ') c_species = c_bro
    WHERE (advt == 'HOBr      ') c_species = c_hobr
    WHERE (advt == 'BrONO2    ') c_species = c_brono2
    WHERE (advt == 'BrCl      ') c_species = c_brcl
    WHERE (advt == 'MeBr      ') c_species = c_mebr
    WHERE (advt == 'HBr       ') c_species = c_hbr
    WHERE (advt == 'CF2ClCFCl2') c_species = c_cf2clcfcl2
    WHERE (advt == 'CHF2Cl    ') c_species = c_chf2cl
    WHERE (advt == 'MeCCl3    ') c_species = c_meccl3
    WHERE (advt == 'CCl4      ') c_species = c_ccl4
    WHERE (advt == 'MeCl      ') c_species = c_mecl
    WHERE (advt == 'CF2ClBr   ') c_species = c_cf2clbr
    WHERE (advt == 'CF3Br     ') c_species = c_cf3br
    WHERE (advt == 'CH2Br2    ') c_species = c_ch2br2
    WHERE (advt == 'C5H8      ') c_species = c_c5h8
    WHERE (advt == 'ISO2      ') c_species = c_iso2
    WHERE (advt == 'ISOOH     ') c_species = c_isooh
    WHERE (advt == 'ISON      ') c_species = c_ison
    WHERE (advt == 'MACR      ') c_species = c_macr
    WHERE (advt == 'MACRO2    ') c_species = c_macro2
    WHERE (advt == 'MACROOH   ') c_species = c_macrooh
    WHERE (advt == 'MPAN      ') c_species = c_mpan
    WHERE (advt == 'HACET     ') c_species = c_hacet
    WHERE (advt == 'MGLY      ') c_species = c_mgly
    WHERE (advt == 'NALD      ') c_species = c_nald
    WHERE (advt == 'HCOOH     ') c_species = c_hcooh
    WHERE (advt == 'MeCO3H    ') c_species = c_meco3h
    WHERE (advt == 'MeCO2H    ') c_species = c_meco2h
    WHERE (advt == 'SO2       ') c_species = c_so2
    WHERE (advt == 'SO3       ') c_species = c_so3
    WHERE (advt == 'DMS       ') c_species = c_dms
    WHERE (advt == 'DMSO      ') c_species = c_dmso
    WHERE (advt == 'Me2S      ') c_species = c_me2s
    WHERE (advt == 'COS       ') c_species = c_cos
    WHERE (advt == 'H2S       ') c_species = c_h2s
    WHERE (advt == 'CS2       ') c_species = c_cs2
    WHERE (advt == 'MSA       ') c_species = c_msa
    WHERE (advt == 'H2SO4     ') c_species = c_h2so4
    WHERE (advt == 'NH3       ') c_species = c_nh3
    WHERE (advt == 'MeOH      ') c_species = c_meoh
    WHERE (advt == 'Monoterp  ') c_species = c_monoterp
    WHERE (advt == 'Sec_Org   ') c_species = c_sec_org
    WHERE (advt == 'SEC_ORG_I ') c_species = c_sec_org_i
    WHERE (advt == 'O(3P)     ') c_species = c_o3p
    WHERE (advt == 'OH        ') c_species = c_oh
    WHERE (advt == 'HO2       ') c_species = c_ho2
    WHERE (advt == 'MeOO      ') c_species = c_meoo
    WHERE (advt == 'EtOO      ') c_species = c_etoo
    WHERE (advt == 'MeCO3     ') c_species = c_meco3
    WHERE (advt == 'n-PrOO    ') c_species = c_proo
    WHERE (advt == 'i-PrOO    ') c_species = c_proo
    WHERE (advt == 'EtCO3     ') c_species = c_etco3
    WHERE (advt == 'MeCOCH2OO ') c_species = c_mecoch2oo
    WHERE (advt == 'O(3P)S    ') c_species = c_o3p
    WHERE (advt == 'O(1D)S    ') c_species = c_o1d
    WHERE (advt == 'OHS       ') c_species = c_oh
    WHERE (advt == 'HO2S      ') c_species = c_ho2
    WHERE (advt == 'N         ') c_species = c_n
    WHERE (advt == 'H         ') c_species = c_h
    WHERE (advt == 'H2        ') c_species = c_h2
    WHERE (advt == 'C4H10     ') c_species = c_c4h10
    WHERE (advt == 'MEK       ') c_species = c_mek
    WHERE (advt == 'C2H4      ') c_species = c_c2h4
    WHERE (advt == 'C3H6      ') c_species = c_c3h6
    WHERE (advt == 'oXYLENE   ') c_species = c_oxylene
    WHERE (advt == 's-BuOOH   ') c_species = c_buooh
    WHERE (advt == 'CH3OH     ') c_species = c_ch3oh
    WHERE (advt == 'MeOH      ') c_species = c_ch3oh
    WHERE (advt == 'GLY       ') c_species = c_gly
    WHERE (advt == 'MEMALD    ') c_species = c_memald
    WHERE (advt == 'MVK       ') c_species = c_mvk
    WHERE (advt == 'MVKOOH    ') c_species = c_mvkooh
    WHERE (advt == 'TOLUENE   ') c_species = c_toluene
    WHERE (advt == 'RNC2H4    ') c_species = c_rnc2h4
    WHERE (advt == 'RNC3H6    ') c_species = c_rnc3h6
    WHERE (advt == 'ORGNIT    ') c_species = c_orgnit
    WHERE (advt == 'NOA       ') c_species = c_noa
    WHERE (advt == 'MeO2NO2   ') c_species = c_meo2no2
    WHERE (advt == 'EtONO2    ') c_species = c_etono2
    WHERE (advt == 'i-PrONO2  ') c_species = c_prono2
    WHERE (advt == 'C2H2      ') c_species = c_c2h2
    WHERE (advt == 'BENZENE   ') c_species = c_benzene
    WHERE (advt == 'APINENE   ') c_species = c_monoterp
    WHERE (advt == 'BPINENE   ') c_species = c_monoterp
    WHERE (advt == 'TBUT2ENE  ') c_species = c_tbut2ene
    WHERE (advt == 'n-PrOH    ') c_species = c_proh
    WHERE (advt == 'i-PrOH    ') c_species = c_proh
    WHERE (advt == 'EtOH      ') c_species = c_etoh
    WHERE (advt == 'EtCO3H    ') c_species = c_etco3h
    WHERE (advt == 'HOCH2CO3  ') c_species = c_hoch2co3
    WHERE (advt == 'HOCH2CO3H ') c_species = c_hoch2co3h
    WHERE (advt == 'HOCH2CH2O2') c_species = c_hoch2ch2o2
    WHERE (advt == 'HOCH2CHO  ') c_species = c_hoch2cho
    WHERE (advt == 'HOC2H4OOH ') c_species = c_hoc2h4ooh
    WHERE (advt == 'HOC2H4NO3 ') c_species = c_hoc2h4no3
    WHERE (advt == 'PHAN      ') c_species = c_phan
    WHERE (advt == 'MeSCH2OO ') c_species = c_ch3sch2oo
    WHERE (advt == 'MeS      ') c_species = c_ch3s
    WHERE (advt == 'MeSO     ') c_species = c_ch3so
    WHERE (advt == 'MeSO2    ') c_species = c_ch3so2
    WHERE (advt == 'MeSO3    ') c_species = c_ch3so3
    WHERE (advt == 'MSIA      ') c_species = c_msia
    WHERE (advt == 'CARB7     ') c_species = c_hacet
    WHERE (advt == 'CARB3     ') c_species = c_gly
    WHERE (advt == 'CARB6     ') c_species = c_mgly
    WHERE (advt == 'UCARB10   ') c_species = c_macr
    WHERE (advt == 'RN10NO3   ') c_species = c_prono2
    WHERE (advt == 'RN13NO3   ') c_species = c_rn13no3
    WHERE (advt == 'RN16NO3   ') c_species = c_rn16no3
    WHERE (advt == 'RN19NO3   ') c_species = c_rn19no3
    WHERE (advt == 'RA13NO3   ') c_species = c_ra13no3
    WHERE (advt == 'RA16NO3   ') c_species = c_ra16no3
    WHERE (advt == 'RA19NO3   ') c_species = c_ra19no3
    WHERE (advt == 'RTX24NO3  ') c_species = c_rtx24no3
    WHERE (advt == 'RN9NO3    ') c_species = c_rn9no3
    WHERE (advt == 'RN12NO3   ') c_species = c_rn12no3
    WHERE (advt == 'RN15NO3   ') c_species = c_rn15no3
    WHERE (advt == 'RN18NO3   ') c_species = c_rn18no3
    WHERE (advt == 'RU14NO3   ') c_species = c_ru14no3
    WHERE (advt == 'RTN28NO3  ') c_species = c_rtn28no3
    WHERE (advt == 'RTN25NO3  ') c_species = c_rtn25no3
    WHERE (advt == 'RTN23NO3  ') c_species = c_rtn23no3
    WHERE (advt == 'RTX28NO3  ') c_species = c_rtx28no3
    WHERE (advt == 'RTX22NO3  ') c_species = c_rtx22no3
    WHERE (advt == 'RN10OOH   ') c_species = c_prooh
    WHERE (advt == 'RN13OOH   ') c_species = c_buooh
    WHERE (advt == 'RN16OOH   ') c_species = c_rn16ooh
    WHERE (advt == 'RN19OOH   ') c_species = c_rn19ooh
    WHERE (advt == 'RN8OOH    ') c_species = c_alkaooh
    WHERE (advt == 'RN11OOH   ') c_species = c_mekooh
    WHERE (advt == 'RN14OOH   ') c_species = c_rn14ooh
    WHERE (advt == 'RN17OOH   ') c_species = c_rn17ooh
    WHERE (advt == 'RU14OOH   ') c_species = c_isooh
    WHERE (advt == 'RU12OOH   ') c_species = c_cri
    WHERE (advt == 'RU10OOH   ') c_species = c_macrooh
    WHERE (advt == 'NRU14OOH  ') c_species = c_nru14ooh
    WHERE (advt == 'NRU12OOH  ') c_species = c_nru12ooh
    WHERE (advt == 'RN9OOH    ') c_species = c_rn9ooh
    WHERE (advt == 'RN12OOH   ') c_species = c_rn12ooh
    WHERE (advt == 'RN15OOH   ') c_species = c_rn15ooh
    WHERE (advt == 'RN18OOH   ') c_species = c_rn18ooh
    WHERE (advt == 'NRN6OOH   ') c_species = c_nrn6ooh
    WHERE (advt == 'NRN9OOH   ') c_species = c_nrn9ooh
    WHERE (advt == 'NRN12OOH  ') c_species = c_nrn12ooh
    WHERE (advt == 'RA13OOH   ') c_species = c_ra13ooh
    WHERE (advt == 'RA16OOH   ') c_species = c_ra16ooh
    WHERE (advt == 'RA19OOH   ') c_species = c_ra19ooh
    WHERE (advt == 'RTN28OOH  ') c_species = c_rtn28ooh
    WHERE (advt == 'NRTN28OOH ') c_species = c_nrtn28ooh
    WHERE (advt == 'RTN26OOH  ') c_species = c_rtn26ooh
    WHERE (advt == 'RTN25OOH  ') c_species = c_rtn25ooh
    WHERE (advt == 'RTN24OOH  ') c_species = c_rtn24ooh
    WHERE (advt == 'RTN23OOH  ') c_species = c_rtn23ooh
    WHERE (advt == 'RTN14OOH  ') c_species = c_rtn14ooh
    WHERE (advt == 'RTN10OOH  ') c_species = c_rtn10ooh
    WHERE (advt == 'RTX28OOH  ') c_species = c_rtx28ooh
    WHERE (advt == 'RTX24OOH  ') c_species = c_rtx24ooh
    WHERE (advt == 'RTX22OOH  ') c_species = c_rtx22ooh
    WHERE (advt == 'NRTX28OOH ') c_species = c_nrtx28ooh
    WHERE (advt == 'CARB14    ') c_species = c_carb14
    WHERE (advt == 'CARB17    ') c_species = c_carb17
    WHERE (advt == 'CARB11A   ') c_species = c_carb11a
    WHERE (advt == 'CARB10    ') c_species = c_carb10
    WHERE (advt == 'CARB13    ') c_species = c_carb13
    WHERE (advt == 'CARB16    ') c_species = c_carb16
    WHERE (advt == 'CARB9     ') c_species = c_carb9
    WHERE (advt == 'CARB12    ') c_species = c_carb12
    WHERE (advt == 'CARB15    ') c_species = c_carb15
    WHERE (advt == 'UCARB12   ') c_species = c_ucarb12
    WHERE (advt == 'NUCARB12  ') c_species = c_nucarb12
    WHERE (advt == 'UDCARB8   ') c_species = c_udcarb8
    WHERE (advt == 'UDCARB11  ') c_species = c_udcarb11
    WHERE (advt == 'UDCARB14  ') c_species = c_udcarb14
    WHERE (advt == 'TNCARB26  ') c_species = c_tncarb26
    WHERE (advt == 'TNCARB10  ') c_species = c_tncarb10
    WHERE (advt == 'TNCARB12  ') c_species = c_tncarb12
    WHERE (advt == 'TNCARB11  ') c_species = c_tncarb11
    WHERE (advt == 'CCARB12   ') c_species = c_ccarb12
    WHERE (advt == 'TNCARB15  ') c_species = c_tncarb15
    WHERE (advt == 'RCOOH25   ') c_species = c_rcooh25
    WHERE (advt == 'TXCARB24  ') c_species = c_txcarb24
    WHERE (advt == 'TXCARB22  ') c_species = c_txcarb22
    WHERE (advt == 'RU12PAN   ') c_species = c_ru12pan
    WHERE (advt == 'RTN26PAN  ') c_species = c_rtn26pan
    WHERE (advt == 'AROH14    ') c_species = c_aroh14
    WHERE (advt == 'AROH17    ') c_species = c_aroh17
    WHERE (advt == 'ARNOH14   ') c_species = c_arnoh14
    WHERE (advt == 'ARNOH17   ') c_species = c_arnoh17
    WHERE (advt == 'ANHY      ') c_species = c_anhy
    WHERE (advt == 'DHPCARB9  ') c_species = c_dhpcarb9
    WHERE (advt == 'HPUCARB12 ') c_species = c_hpucarb12
    WHERE (advt == 'HUCARB9   ') c_species = c_hucarb9
    WHERE (advt == 'IEPOX     ') c_species = c_iepox
    WHERE (advt == 'HMML      ') c_species = c_hmml
    WHERE (advt == 'DHPR12OOH ') c_species = c_dhpr12ooh
    WHERE (advt == 'DHCARB9   ') c_species = c_dhcarb9
    WHERE (advt == 'RU12NO3   ') c_species = c_ru12no3
    WHERE (advt == 'RU10NO3   ') c_species = c_ru10no3
    WHERE (advt == 'RAROH14   ') c_species = c_cri
    WHERE (advt == 'RAROH17   ') c_species = c_cri
    WHERE (advt == 'PASSIVE O3') c_species = 1.0
    WHERE (advt == 'AGE OF AIR') c_species = 1.0
    !$acc kernels
    !$acc loop gang vector independent
    do idx_1 = LBOUND(c_na_species, dim=1), UBOUND(c_na_species, dim=1), 1
      c_na_species(idx_1) = 0.0
    enddo
    !$acc end kernels
    ! PSyclone CodeBlock (unsupported code) reason:
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    !  - PSyclone doesn't yet support reference to imported symbols inside WHERE clauses.
    WHERE (nadvt == 'OH        ') c_na_species = c_oh
    WHERE (nadvt == 'HO2       ') c_na_species = c_ho2
    WHERE (nadvt == 'O(3P)     ') c_na_species = c_o3p
    WHERE (nadvt == 'O3P       ') c_na_species = c_o3p
    WHERE (nadvt == 'O(1D)     ') c_na_species = c_o1d
    WHERE (nadvt == 'O1D       ') c_na_species = c_o1d
    WHERE (nadvt == 'MeOO      ') c_na_species = c_meoo
    WHERE (nadvt == 'EtOO      ') c_na_species = c_etoo
    WHERE (nadvt == 'MeCO3     ') c_na_species = c_meco3
    WHERE (nadvt == 'n-PrOO    ') c_na_species = c_proo
    WHERE (nadvt == 'i-PrOO    ') c_na_species = c_proo
    WHERE (nadvt == 'EtCO3     ') c_na_species = c_etco3
    WHERE (nadvt == 'MeCOCH2OO ') c_na_species = c_mecoch2oo
    WHERE (nadvt == 'O(3P)S    ') c_na_species = c_o3p
    WHERE (nadvt == 'O(1D)S    ') c_na_species = c_o1d
    WHERE (nadvt == 'OHS       ') c_na_species = c_oh
    WHERE (nadvt == 'HO2S      ') c_na_species = c_ho2
    WHERE (nadvt == 's-BuOO    ') c_na_species = c_buoo
    WHERE (nadvt == 'MEKO2     ') c_na_species = c_meko2
    WHERE (nadvt == 'HOC2H4O2  ') c_na_species = c_hoc2h4o2
    WHERE (nadvt == 'HOC3H6O2  ') c_na_species = c_hoc3h6o2
    WHERE (nadvt == 'OXYL1     ') c_na_species = c_oxyl1
    WHERE (nadvt == 'MEMALD1   ') c_na_species = c_memald1
    WHERE (nadvt == 'HOIPO2    ') c_na_species = c_hoipo2
    WHERE (nadvt == 'HOMVKO2   ') c_na_species = c_homvko2
    WHERE (nadvt == 'TOLP1     ') c_na_species = c_tolp1
    WHERE (nadvt == 'ISO2      ') c_na_species = c_iso2
    WHERE (nadvt == 'MACRO2    ') c_na_species = c_macro2
    WHERE (nadvt == 'HOCH2CO3  ') c_na_species = c_hoch2co3
    WHERE (nadvt == 'HOCH2CH2O2') c_na_species = c_hoch2ch2o2
    WHERE (nadvt == 'RN10O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN13O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN16O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN19O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN13AO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RN16AO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RA13O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RA16O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RA19AO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RA19CO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RN9O2     ') c_na_species = c_cri
    WHERE (nadvt == 'RN12O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN15O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN18O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN15AO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RN18AO2   ') c_na_species = c_cri
    WHERE (nadvt == 'RN8O2     ') c_na_species = c_cri
    WHERE (nadvt == 'RN11O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN14O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RN17O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RU14O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RU12O2    ') c_na_species = c_cri
    WHERE (nadvt == 'RU10O2    ') c_na_species = c_cri
    WHERE (nadvt == 'NRN6O2    ') c_na_species = c_cri
    WHERE (nadvt == 'NRN9O2    ') c_na_species = c_cri
    WHERE (nadvt == 'NRN12O2   ') c_na_species = c_cri
    WHERE (nadvt == 'NRU14O2   ') c_na_species = c_cri
    WHERE (nadvt == 'NRU12O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN28O2   ') c_na_species = c_cri
    WHERE (nadvt == 'NRTN28O2  ') c_na_species = c_cri
    WHERE (nadvt == 'RTN26O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN25O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN24O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN23O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN14O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTN10O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTX28O2   ') c_na_species = c_cri
    WHERE (nadvt == 'NRTX28O2  ') c_na_species = c_cri
    WHERE (nadvt == 'RTX24O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RTX22O2   ') c_na_species = c_cri
    WHERE (nadvt == 'RO2       ') c_na_species = c_cri
    WHERE (nadvt == 'DHPR12O2  ') c_na_species = c_dhpr12o2
    WHERE (nadvt == 'RU10AO2   ') c_na_species = c_ru10ao2
    WHERE (nadvt == 'MACO3     ') c_na_species = c_maco3
    n_h2o = 0
    n_cl = 0
    n_sx = 0
    n_n2o5 = 0
    n_ox = 0
    n_o3 = 0
    n_o3s = 0
    n_nox = 0
    n_no = 0
    n_no2 = 0
    n_no3 = 0
    n_ho2no2 = 0
    n_hono2 = 0
    n_pan = 0
    n_ison = 0
    n_orgnit = 0
    n_rnc2h4 = 0
    n_rnc3h6 = 0
    n_h2o2 = 0
    n_ch4 = 0
    n_h2 = 0
    n_dms = 0
    n_so2 = 0
    n_so3 = 0
    n_h2so4 = 0
    n_sec_org = 0
    n_sec_org_i = 0
    n_passive = 0
    n_hcho = 0
    n_c2h6 = 0
    n_mecho = 0
    n_c3h8 = 0
    n_me2co = 0
    n_c5h8 = 0
    n_mgly = 0
    n_c4h10 = 0
    n_mek = 0
    n_c2h4 = 0
    n_c3h6 = 0
    n_oxylene = 0
    n_ch3oh = 0
    n_gly = 0
    n_mvk = 0
    n_toluene = 0
    n_co = 0
    n_n2o = 0
    n_cf2cl2 = 0
    n_cfcl3 = 0
    n_mebr = 0
    n_nh3 = 0
    do m = 1, jpctr, 1
      if (advt(m) == 'Ox        ') then
        n_ox = m
        n_ox1 = m
      else
        if (advt(m) == 'O3        ') then
          n_o3 = m
        else
          if (advt(m) == 'O3S       ') then
            n_o3s = m
          else
            if (advt(m) == 'Sx        ') then
              n_sx = m
              n_sx1 = m
            else
              if (advt(m) == 'NOx       ') then
                n_nox = m
              else
                if (advt(m) == 'NO        ') then
                  n_no = m
                else
                  if (advt(m) == 'NO2       ') then
                    n_no2 = m
                  else
                    if (advt(m) == 'NO3       ') then
                      n_no3 = m
                    else
                      if (advt(m) == 'N2O5      ') then
                        n_n2o5 = m
                      else
                        if (advt(m) == 'HO2NO2    ') then
                          n_ho2no2 = m
                        else
                          if (advt(m) == 'HONO2     ') then
                            n_hono2 = m
                          else
                            if (advt(m) == 'PAN       ') then
                              n_pan = m
                            else
                              if (advt(m) == 'ISON      ') then
                                n_ison = m
                              else
                                if (advt(m) == 'ORGNIT    ') then
                                  n_orgnit = m
                                else
                                  if (advt(m) == 'RNC2H4    ') then
                                    n_rnc2h4 = m
                                  else
                                    if (advt(m) == 'RNC3H6    ') then
                                      n_rnc3h6 = m
                                    else
                                      if (advt(m) == 'H2O2      ') then
                                        n_h2o2 = m
                                      else
                                        if (advt(m) == 'CH4       ') then
                                          n_ch4 = m
                                        else
                                          if (advt(m) == 'H2O       ') then
                                            n_h2o = m
                                          else
                                            if (advt(m) == 'Cl        ') then
                                              n_cl = m
                                            else
                                              if (advt(m) == 'H2        ') then
                                                n_h2 = m
                                              else
                                                if (advt(m) == 'DMS       ') then
                                                  n_dms = m
                                                else
                                                  if (advt(m) == 'SO2       ') then
                                                    n_so2 = m
                                                  else
                                                    if (advt(m) == 'SO3       ') then
                                                      n_so3 = m
                                                    else
                                                      if (advt(m) == 'H2SO4     ') then
                                                        n_h2so4 = m
                                                      else
                                                        if (advt(m) == 'Sec_Org   ') then
                                                          n_sec_org = m
                                                        else
                                                          if (advt(m) == 'SEC_ORG_I ') then
                                                            n_sec_org_i = m
                                                          else
                                                            if (advt(m) == 'BrO       ') then
                                                              n_bro = m
                                                            else
                                                              if (advt(m) == 'HCl       ') then
                                                                n_hcl = m
                                                              else
                                                                if (advt(m) == 'CFCl3     ') then
                                                                  n_cfcl3 = m
                                                                else
                                                                  if (advt(m) == 'CF2Cl2    ') then
                                                                    n_cf2cl2 = m
                                                                  else
                                                                    if (advt(m) == 'O(1D)     ') then
                                                                      n_o1d = m
                                                                    else
                                                                      if (advt(m) == 'HCHO      ') then
                                                                        n_hcho = m
                                                                      else
                                                                        if (advt(m) == 'C2H6      ') then
                                                                          n_c2h6 = m
                                                                        else
                                                                          if (advt(m) == 'MeCHO     ') then
                                                                            n_mecho = m
                                                                          else
                                                                            if (advt(m) == 'C3H8      ') then
                                                                              n_c3h8 = m
                                                                            else
                                                                              if (advt(m) == 'Me2CO     ') then
                                                                                n_me2co = m
                                                                              else
                                                                                if (advt(m) == 'C5H8      ') then
                                                                                  n_c5h8 = m
                                                                                else
                                                                                  if (advt(m) == 'MGLY      ') then
                                                                                    n_mgly = m
                                                                                  else
                                                                                    if (advt(m) == 'C4H10     ') then
                                                                                      n_c4h10 = m
                                                                                    else
                                                                                      if (advt(m) == 'MEK       ') then
                                                                                        n_mek = m
                                                                                      else
                                                                                        if (advt(m) == 'C2H4      ') then
                                                                                          n_c2h4 = m
                                                                                        else
                                                                                          if (advt(m) == 'C3H6      ') then
                                                                                            n_c3h6 = m
                                                                                          else
                                                                                            if (advt(m) == 'oXYLENE   ') then
                                                                                              n_oxylene = m
                                                                                            else
                                                                                              if (advt(m) == 'CH3OH     ') then
                                                                                                n_ch3oh = m
                                                                                              else
                                                                                                if (advt(m) == 'GLY       ') then
                                                                                                  n_gly = m
                                                                                                else
                                                                                                  if (advt(m) == 'MVK       ') then
                                                                                                    n_mvk = m
                                                                                                  else
                                                                                                    if (advt(m) == 'TOLUENE   ') then
                                                                                                      n_toluene = m
                                                                                                    else
                                                                                                      if (advt(m) == 'CO        ') then
                                                                                                        n_co = m
                                                                                                      else
                                                                                                        if (advt(m) == 'MeBr      ') then
                                                                                                          n_mebr = m
                                                                                                        else
                                                                                                          if (advt(m) == 'N2O       ') then
                                                                                                            n_n2o = m
                                                                                                          else
                                                                                                            if (advt(m) == 'NH3       ') then
                                                                                                              n_nh3 = m
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
    do i = 1, jpctr, 1
      if (advt(i) == 'Ox        ' .OR. advt(i) == 'O3        ') then
        n_o3 = i
      end if
    enddo
    if (n_o3 == 0 .AND. (.NOT.(ukca_config%l_ukca_offline .OR. ukca_config%l_ukca_offline_be))) then
      cmessage(:) = 'Ozone not found among chemical tracers.'
      errcode = 1
      call ereport('UKCA_CSPECIES', errcode, cmessage)
    end if
    !$acc kernels
    nn_o3 = 0
    nn_o3s = 0
    nn_oh = 0
    nn_ho2 = 0
    nn_h2o2 = 0
    nn_no = 0
    nn_no2 = 0
    nn_o1d = 0
    nn_o3p = 0
    nn_meoo = 0
    nn_meco3 = 0
    nn_etoo = 0
    nn_etco3 = 0
    nn_nproo = 0
    nn_nprooh = 0
    nn_iproo = 0
    nn_iprooh = 0
    nn_mecoch2oo = 0
    nn_so2 = 0
    nn_so3 = 0
    nn_h2so4 = 0
    nn_ohs = 0
    nn_ho2s = 0
    nn_o1ds = 0
    nn_o3ps = 0
    nn_n = 0
    nn_meo = 0
    nn_cl = 0
    nn_clo = 0
    nn_cl2o2 = 0
    nn_hcl = 0
    nn_br = 0
    nn_bro = 0
    nn_buoo = 0
    nn_meko2 = 0
    nn_hoc2h4o2 = 0
    nn_hoc3h6o2 = 0
    nn_oxyl1 = 0
    nn_memald1 = 0
    nn_hoipo2 = 0
    nn_homvko2 = 0
    nn_tolp1 = 0
    nn_nh3 = 0
    !$acc loop seq
    do m = 1, jpspec, 1
      if (speci(m) == 'O3        ') then
        nn_o3 = m
      else
        if (speci(m) == 'O3S       ') then
          nn_o3s = m
        else
          if (speci(m) == 'OH        ') then
            nn_oh = m
          else
            if (speci(m) == 'HO2       ') then
              nn_ho2 = m
            else
              if (speci(m) == 'H2O2      ') then
                nn_h2o2 = m
              else
                if (speci(m) == 'NO        ') then
                  nn_no = m
                else
                  if (speci(m) == 'NO2       ') then
                    nn_no2 = m
                  else
                    if (speci(m) == 'O(1D)     ') then
                      nn_o1d = m
                    else
                      if (speci(m) == 'O(3P)     ') then
                        nn_o3p = m
                      else
                        if (speci(m) == 'MeOO      ') then
                          nn_meoo = m
                        else
                          if (speci(m) == 'MeCO3     ') then
                            nn_meco3 = m
                          else
                            if (speci(m) == 'EtOO      ') then
                              nn_etoo = m
                            else
                              if (speci(m) == 'EtCO3     ') then
                                nn_etco3 = m
                              else
                                if (speci(m) == 'n-PrOO    ') then
                                  nn_nproo = m
                                else
                                  if (speci(m) == 'n-PrOOH   ') then
                                    nn_nprooh = m
                                  else
                                    if (speci(m) == 'i-PrOO    ') then
                                      nn_iproo = m
                                    else
                                      if (speci(m) == 'i-PrOOH   ') then
                                        nn_iprooh = m
                                      else
                                        if (speci(m) == 'MeCOCH2OO ') then
                                          nn_mecoch2oo = m
                                        else
                                          if (speci(m) == 'NO3       ') then
                                            nn_no3 = m
                                          else
                                            if (speci(m) == 'N2O5      ') then
                                              nn_n2o5 = m
                                            else
                                              if (speci(m) == 'HO2NO2    ') then
                                                nn_ho2no2 = m
                                              else
                                                if (speci(m) == 'HONO2     ') then
                                                  nn_hono2 = m
                                                else
                                                  if (speci(m) == 'CH4       ') then
                                                    nn_ch4 = m
                                                  else
                                                    if (speci(m) == 'SO2       ') then
                                                      nn_so2 = m
                                                    else
                                                      if (speci(m) == 'SO3       ') then
                                                        nn_so3 = m
                                                      else
                                                        if (speci(m) == 'H2SO4     ') then
                                                          nn_h2so4 = m
                                                        else
                                                          if (speci(m) == 'O(1D)S    ') then
                                                            nn_o1ds = m
                                                          else
                                                            if (speci(m) == 'O(3P)S    ') then
                                                              nn_o3ps = m
                                                            else
                                                              if (speci(m) == 'OHS       ') then
                                                                nn_ohs = m
                                                              else
                                                                if (speci(m) == 'HO2S      ') then
                                                                  nn_ho2s = m
                                                                else
                                                                  if (speci(m) == 'N         ') then
                                                                    nn_n = m
                                                                  else
                                                                    if (speci(m) == 'MeO       ') then
                                                                      nn_meo = m
                                                                    else
                                                                      if (speci(m) == 'Cl        ') then
                                                                        nn_cl = m
                                                                      else
                                                                        if (speci(m) == 'ClO       ') then
                                                                          nn_clo = m
                                                                        else
                                                                          if (speci(m) == 'Cl2O2     ') then
                                                                            nn_cl2o2 = m
                                                                          else
                                                                            if (speci(m) == 'HCl       ') then
                                                                              nn_hcl = m
                                                                            else
                                                                              if (speci(m) == 'Br        ') then
                                                                                nn_br = m
                                                                              else
                                                                                if (speci(m) == 'BrO       ') then
                                                                                  nn_bro = m
                                                                                else
                                                                                  if (speci(m) == 's-BuOO    ') then
                                                                                    nn_buoo = m
                                                                                  else
                                                                                    if (speci(m) == 'MEKO2     ') then
                                                                                      nn_meko2 = m
                                                                                    else
                                                                                      if (speci(m) == 'HOC2H4O2  ') then
                                                                                        nn_hoc2h4o2 = m
                                                                                      else
                                                                                        if (speci(m) == 'HOC3H6O2  ') then
                                                                                          nn_hoc3h6o2 = m
                                                                                        else
                                                                                          if (speci(m) == 'OXYL1     ') then
                                                                                            nn_oxyl1 = m
                                                                                          else
                                                                                            if (speci(m) == 'MEMALD1   ') then
                                                                                              nn_memald1 = m
                                                                                            else
                                                                                              if (speci(m) == 'HOIPO2    ') then
                                                                                                nn_hoipo2 = m
                                                                                              else
                                                                                                if (speci(m) == 'HOMVKO2   ') then
                                                                                                  nn_homvko2 = m
                                                                                                else
                                                                                                  if (speci(m) == 'TOLP1     ') then
                                                                                                    nn_tolp1 = m
                                                                                                  else
                                                                                                    if (speci(m) == 'CFCl3     ') then
                                                                                                      nn_cfcl3 = m
                                                                                                    else
                                                                                                      if (speci(m) == 'CF2Cl2    ') then
                                                                                                        nn_cf2cl2 = m
                                                                                                      else
                                                                                                        if (speci(m) == 'NH3       ') then
                                                                                                          nn_nh3 = m
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
    if (printstatus > prstatus_oper) then
      call umprint('Species indices:', src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'O3 = ', nn_o3, ' OH = ', nn_oh
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'HO2 = ', nn_ho2, ' NO = ', nn_no
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'NO2 = ', nn_no2, ' O(1D) = ', nn_o1d
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'MeOO = ', nn_meoo, ' MeCO3 = ', nn_meco3
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'EtOO = ', nn_etoo, ' EtCO3 = ', nn_etco3
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'n-Proo = ', nn_nproo, 'n-Prooh = ', nn_nprooh
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'i-Proo = ', nn_iproo, 'i-Prooh = ', nn_iprooh
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'MeCOCH2OO = ', nn_mecoch2oo
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'SO2: ', nn_so2, ' H2O2: ', nn_h2o2
      call umprint(ummessage, src='ukca_cspecies')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Write_Stmt
      WRITE(umMessage, '(A,1X,I6,A,1X,I6)') 'SO3: ', nn_so3, ' H2SO4: ', nn_h2so4
      call umprint(ummessage, src='ukca_cspecies')
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_calc_cspecies

end module ukca_cspecies
