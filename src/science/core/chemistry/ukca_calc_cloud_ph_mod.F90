module ukca_calc_cloud_ph_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CALC_CLOUD_PH_MOD'
  private

  public :: ukca_calc_cloud_ph

  contains
  subroutine ukca_calc_cloud_ph(tracer, row_length, rows, model_levels, pres, temp, h_plus_3d_arr)
    use ukca_constants, only : h_plus, m_so2
    use ukca_config_constants_mod, only : avogadro, boltzmann
    use ukca_cspecies, only : c_species, n_so2
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    real, parameter :: kgpcm3_to_ugpm3 = 1.0e15
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    real, dimension(row_length,rows,model_levels), intent(in) :: pres
    real, dimension(row_length,rows,model_levels), intent(in) :: temp
    real, dimension(:,:,:,:), intent(in) :: tracer
    real, dimension(:,:,:), intent(inout) :: h_plus_3d_arr
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    real :: so2_conc
    real :: so2_conc_micro
    real :: so2_conc_micro_log
    real :: ph_val
    real :: cur_h_plus
    real :: new_h_plus
    real, dimension(row_length,rows,model_levels) :: aird
    real :: fit_1
    real :: fit_2
    real :: fit_3
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CALC_CLOUD_PH'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    !$acc kernels
    fit_1 = ukca_config%ph_fit_coeff_a
    fit_2 = ukca_config%ph_fit_coeff_b
    fit_3 = ukca_config%ph_fit_intercept
    aird(:,:,:) = pres(:,:,:) / (temp(:,:,:) * boltzmann * 1.0e6)
    !$acc loop seq
    do k = 1, model_levels, 1
      !$acc loop seq
      do j = 1, rows, 1
        !$acc loop seq
        do i = 1, row_length, 1
          so2_conc = tracer(i,j,k,n_so2) / c_species(n_so2)
          so2_conc_micro = so2_conc * aird(i,j,k) * m_so2 * 1.0e-3 * kgpcm3_to_ugpm3 / avogadro
          so2_conc_micro_log = LOG10(so2_conc_micro)
          ph_val = fit_1 * so2_conc_micro_log * so2_conc_micro_log + fit_2 * so2_conc_micro_log + fit_3
          cur_h_plus = 10.0 ** (-1.0 * ph_val)
          if (ph_val < 2.0) then
            ph_val = 2.00
          end if
          if (ph_val > 8.3) then
            ph_val = 8.30
          end if
          new_h_plus = 10.0 ** (-1.0 * ph_val)
          h_plus_3d_arr(i,j,k) = new_h_plus
        enddo
      enddo
    enddo
    !$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if

  end subroutine ukca_calc_cloud_ph

end module ukca_calc_cloud_ph_mod
