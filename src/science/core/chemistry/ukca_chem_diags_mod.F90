module ukca_chem_diags_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEM_DIAGS_MOD'
  public

  contains
  subroutine ukca_chem_diags(row_length, rows, model_levels, nat_psc, trop_ch4_mol, trop_o3_mol, trop_oh_mol, strat_ch4_mol, strat_ch4loss, atm_ch4_mol, atm_co_mol, atm_n2o_mol, atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol, atm_h2_mol, so4_sa, h_plus_3d_arr, len_stashwork, stashwork)
    use ukca_um_legacy_mod, only : copydiag_3d, len_stlist, num_stash_levels, sf, si, si_last, stash_levels, stashcode_ukca_atmos_cfc11, stashcode_ukca_atmos_cfc12, stashcode_ukca_atmos_ch3br, stashcode_ukca_atmos_ch4, stashcode_ukca_atmos_co, stashcode_ukca_atmos_h2, stashcode_ukca_atmos_n2o, stashcode_ukca_chem_diag, stashcode_ukca_h_plus, stashcode_ukca_nat, stashcode_ukca_so4_sad, stashcode_ukca_strat_ch4, stashcode_ukca_strt_ch4_lss, stashcode_ukca_trop_ch4, stashcode_ukca_trop_o3, stashcode_ukca_trop_oh, stindex, stlist, ukca_diag_sect
    use ereport_mod, only : ereport
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    real, dimension(row_length,rows,model_levels), intent(in) :: nat_psc
    real, dimension(row_length,rows,model_levels), intent(in) :: trop_ch4_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: trop_o3_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: trop_oh_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: strat_ch4_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: strat_ch4loss
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_ch4_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_co_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_n2o_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_cf2cl2_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_cfcl3_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_mebr_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: atm_h2_mol
    real, dimension(row_length,rows,model_levels), intent(in) :: so4_sa
    real, dimension(row_length,rows,model_levels), intent(in) :: h_plus_3d_arr
    integer, intent(in) :: len_stashwork
    real, dimension(len_stashwork), intent(inout) :: stashwork
    integer :: item
    integer :: icode
    integer :: im_index
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEM_DIAGS'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    icode = 0
    im_index = 1
    item = stashcode_ukca_nat - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), nat_psc(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_trop_ch4 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), trop_ch4_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_trop_o3 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), trop_o3_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_trop_oh - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), trop_oh_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_strat_ch4 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), strat_ch4_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_strt_ch4_lss - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), strat_ch4loss(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_ch4 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_ch4_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_co - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_co_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_n2o - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_n2o_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_cfc12 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_cf2cl2_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_cfc11 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_cfcl3_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_ch3br - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_mebr_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_atmos_h2 - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), atm_h2_mol(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_so4_sad - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), so4_sa(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    item = stashcode_ukca_h_plus - 1000 * stashcode_ukca_chem_diag
    if (sf(item,ukca_diag_sect)) then
      call copydiag_3d(stashwork(si(item,ukca_diag_sect,im_index):si_last(item,ukca_diag_sect,im_index)), h_plus_3d_arr(:,:,:), row_length, rows, model_levels, stlist(:,stindex(1,item,ukca_diag_sect,im_index)), len_stlist, stash_levels, num_stash_levels + 1)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    return

  end subroutine ukca_chem_diags

end module ukca_chem_diags_mod
