module ukca_chemistry_cleanup_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CHEMISTRY_CLEANUP_MOD'
  public

  contains
  subroutine ukca_chemistry_cleanup(row_length, rows, model_levels, ntracers, secs_per_step, pres, drain, crain, latitude, env_ozone3d, qcf, r_theta_levels, mass, z_top_of_model, shno3_3d, nat_psc, tracer)
    use asad_mod, only : jpctr
    use ukca_config_specification_mod, only : i_top_1lev, i_top_2levh2o, i_top_bc, ukca_config
    use ukca_conserve_mod, only : ukca_conserve
    use ukca_cspecies, only : n_h2o, n_hono2, nn_cl
    use ukca_sediment_mod, only : ukca_sediment
    use ukca_stratf_mod, only : ukca_stratf
    use ukca_topboundary_mod, only : ukca_top_boundary
    use ukca_tropopause, only : l_stratosphere
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    logical, parameter :: after_chem = .false.
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: ntracers
    real, intent(in) :: secs_per_step
    real, dimension(row_length,rows,model_levels), intent(in) :: pres
    real, dimension(row_length,rows,model_levels), intent(in) :: drain
    real, dimension(row_length,rows,model_levels), intent(in) :: crain
    real, dimension(row_length,rows), intent(in) :: latitude
    real, dimension(row_length,rows,model_levels), intent(in) :: env_ozone3d
    real, dimension(row_length,rows,model_levels), intent(in) :: qcf
    real, dimension(row_length,rows,0:model_levels), intent(in) :: r_theta_levels
    real, dimension(row_length,rows,model_levels), intent(in) :: mass
    real, intent(in) :: z_top_of_model
    real, dimension(row_length,rows,model_levels), intent(inout) :: shno3_3d
    real, dimension(row_length,rows,model_levels), intent(out) :: nat_psc
    real, dimension(row_length,rows,model_levels,ntracers), intent(inout) :: tracer
    integer :: i
    integer :: j
    integer :: k
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CHEMISTRY_CLEANUP'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (ukca_config%l_ukca_offline_be) then
      return
    end if
    if (nn_cl > 0) then
      call ukca_conserve(row_length, rows, model_levels, ntracers, tracer, pres, drain, crain, after_chem)
    end if
    if (ukca_config%l_ukca_strat .OR. ukca_config%l_ukca_stratcfc .OR. ukca_config%l_ukca_strattrop .OR. ukca_config%l_ukca_cristrat) then
      if (ukca_config%l_ukca_het_psc) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(3)
        do idx = LBOUND(tracer, dim=3), UBOUND(tracer, dim=3), 1
          do idx_1 = LBOUND(tracer, dim=2), UBOUND(tracer, dim=2), 1
            do idx_2 = LBOUND(tracer, dim=1), UBOUND(tracer, dim=1), 1
              tracer(idx_2,idx_1,idx,n_hono2) = tracer(idx_2,idx_1,idx,n_hono2) - shno3_3d(idx_2,idx_1,idx)
            enddo
          enddo
        enddo
        !NO!$acc end kernels
        call ukca_sediment(rows, row_length, model_levels, shno3_3d, qcf, r_theta_levels, mass, secs_per_step, l_stratosphere(:,:,:))
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(3)
        do idx_3 = LBOUND(tracer, dim=3), UBOUND(tracer, dim=3), 1
          do idx_4 = LBOUND(tracer, dim=2), UBOUND(tracer, dim=2), 1
            do idx_5 = LBOUND(tracer, dim=1), UBOUND(tracer, dim=1), 1
              tracer(idx_5,idx_4,idx_3,n_hono2) = tracer(idx_5,idx_4,idx_3,n_hono2) + shno3_3d(idx_5,idx_4,idx_3)
            enddo
          enddo
        enddo
        !NO!$acc end kernels
      end if
      if (ukca_config%i_ukca_topboundary == i_top_2levh2o) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(3)
        do i = 1, rows, 1
          do j = 1, row_length, 1
            do k = 1, ntracers, 1
              if (k /= n_h2o) then
                tracer(j,i,model_levels,k) = tracer(j,i,model_levels - 2,k)
                tracer(j,i,model_levels - 1,k) = tracer(j,i,model_levels - 2,k)
              end if
            enddo
          enddo
        enddo
        !NO!$acc end kernels
      else
        if (ukca_config%i_ukca_topboundary == i_top_1lev) then
          !NO!$acc kernels
          !NO!$acc loop gang vector independent collapse(3)
          do i = 1, rows, 1
            do j = 1, row_length, 1
              do k = 1, ntracers, 1
                tracer(j,i,model_levels,k) = tracer(j,i,model_levels - 1,k)
              enddo
            enddo
          enddo
          !NO!$acc end kernels
        else
          if (ukca_config%i_ukca_topboundary >= i_top_bc) then
            if (z_top_of_model > 85500.0 .OR. z_top_of_model < 79000.0) then
              errcode = 25
              cmessage(:) = 'Can only impose a top boundary condition at 85 km.'
              call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
            end if
            if (.NOT.ALL(ABS(latitude(row_length,:) - latitude(1,:)) < EPSILON(0.0))) then
              errcode = 1
              cmessage(:) = 'Can only impose top boundary condition if latitude is constant on rows'
              call ereport(ModuleName // ':' // RoutineName, errcode, cmessage)
            end if
            call ukca_top_boundary(row_length, rows, model_levels, ntracers, latitude(1,:), tracer)
          end if
        end if
      end if
      if (ukca_config%l_ukca_het_psc) then
        !NO!$acc kernels
        !NO!$acc loop gang vector independent collapse(3)
        do idx_6 = LBOUND(nat_psc, dim=3), UBOUND(nat_psc, dim=3), 1
          do idx_7 = LBOUND(nat_psc, dim=2), UBOUND(nat_psc, dim=2), 1
            do idx_8 = LBOUND(nat_psc, dim=1), UBOUND(nat_psc, dim=1), 1
              nat_psc(idx_8,idx_7,idx_6) = shno3_3d(idx_8,idx_7,idx_6)
            enddo
          enddo
        enddo
        !NO!$acc end kernels
      end if
    else
      if (.NOT.ukca_config%l_ukca_offline) then
        call ukca_stratf(row_length, rows, model_levels, jpctr, env_ozone3d, tracer(:,:,:,:jpctr))
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine ukca_chemistry_cleanup

end module ukca_chemistry_cleanup_mod
