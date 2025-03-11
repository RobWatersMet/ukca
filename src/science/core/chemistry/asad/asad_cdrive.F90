module asad_cdrive_mod
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ASAD_CDRIVE_MOD'
  public

  contains
  subroutine asad_cdrive(cdot, ftr, pp, pt, pq, co2_1d, cld_f, cld_l, ix, jy, nlev, dryrt, wetrt, rc_het, prt, n_points, have_nat, stratflag, h_plus_1d_arr)
    use asad_mod, only : cdt, co2, ctype, f, fdot, jpcspf, jpdd, jpdw, jpif, jppj, jpspec, jsubs, linfam, lvmr, madvtr, method, moffam, ncsteps, ndepd, ndepw, nfphot, nit0, nitfg, nodd, p, t, tnd, wp
    use ukca_hetero_mod, only : ukca_hetero, ukca_solidphase
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use nvtx
    use ereport_mod, only : ereport
    use umPrintMgr, only : ummessage, umprint
    use errormessagelength_mod, only : errormessagelength
    use asad_bedriv_mod, only : asad_bedriv
    use asad_bimol_mod, only : asad_bimol
    use asad_diffun_mod, only : asad_diffun
    use asad_hetero_mod, only : asad_hetero
    use asad_impact_mod, only : asad_impact
    use asad_jac_mod, only : asad_jac
    use asad_spmjpdriv_mod, only : asad_spmjpdriv
    use asad_totnud_mod, only : asad_totnud
    use asad_trimol_mod, only : asad_trimol
    use ukca_drydep_mod, only : ukca_drydep
    use ukca_wetdep_mod, only : ukca_wetdep
    use ukca_photol_mod, only : ukca_photol
    use asad_posthet_mod, only : asad_posthet
    use asad_ftoy_mod, only : asad_ftoy
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    integer, intent(in) :: ix
    integer, intent(in) :: jy
    integer, intent(in) :: nlev
    real, dimension(n_points,jppj), intent(in) :: prt
    real, dimension(n_points,jpdd), intent(in) :: dryrt
    real, dimension(n_points,jpdw), intent(in) :: wetrt
    real, dimension(n_points,2), intent(in) :: rc_het
    real, dimension(n_points), intent(in) :: pp
    real, dimension(n_points), intent(in) :: pt
    real, dimension(n_points), intent(in) :: pq
    real, dimension(n_points), intent(in) :: co2_1d
    real, dimension(n_points), intent(in) :: cld_f
    real, dimension(n_points), intent(in) :: cld_l
    logical, dimension(n_points), intent(in) :: have_nat
    logical, dimension(n_points), intent(in) :: stratflag
    real, dimension(n_points), intent(in) :: h_plus_1d_arr
    real, dimension(n_points,jpcspf), intent(inout) :: ftr
    real, dimension(n_points,jpcspf), intent(out) :: cdot
    integer :: errcode
    integer :: jtr
    integer :: jl
    integer :: js
    integer :: nl
    integer :: ifam
    integer :: itr
    integer :: iodd
    integer :: num_iter
    logical :: gfirst
    logical :: gphot
    logical, save :: first_call = .true.
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_CDRIVE'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(2)
    do jtr = 1, jpcspf, 1
      do jl = 1, n_points, 1
        cdot(jl,jtr) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do jl = 1, n_points, 1
      p(jl) = pp(jl)
      t(jl) = pt(jl)
    enddo
    !NO!$acc loop gang vector independent
    do jl = 1, n_points, 1
      wp(jl) = pq(jl)
      co2(jl) = co2_1d(jl)
    enddo
    !NO!$acc end kernels
    call asad_totnud(n_points)
    if (lvmr) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do jtr = 1, jpcspf, 1
        do jl = 1, n_points, 1
          ftr(jl,jtr) = ftr(jl,jtr) * tnd(jl)
          f(jl,jtr) = ftr(jl,jtr)
        enddo
      enddo
      !NO!$acc end kernels
    else
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do jtr = 1, jpcspf, 1
        do jl = 1, n_points, 1
          f(jl,jtr) = ftr(jl,jtr)
        enddo
      enddo
      !NO!$acc end kernels
    end if
    call asad_bimol(n_points)
    call asad_trimol(n_points)
    if (ukca_config%l_ukca_nr_aqchem .OR. ukca_config%l_ukca_trophet) then
      call asad_hetero(n_points, cld_f, cld_l, rc_het, h_plus_1d_arr)
    end if
    if (ndepw /= 0) then
      call ukca_wetdep(wetrt, n_points)
    end if
    if (ndepd /= 0) then
      call ukca_drydep(nlev, dryrt, n_points)
    end if
    gphot = .true.
    if (method /= 0) then
      do js = 1, ncsteps, 1
        jsubs = js
        gfirst = jsubs == 1
        gphot = gfirst
        if (nfphot /= 0 .AND. (.NOT.gfirst)) then
          gphot = MOD(jsubs - 1, nfphot) == 0
        end if
        nl = n_points
        num_iter = 0
        if (method == 1) then
          if (gphot) then
            call ukca_photol(prt, n_points)
          end if
          if (ukca_config%l_ukca_het_psc) then
            call ukca_hetero(n_points, have_nat, stratflag)
          end if
          call asad_ftoy(gfirst, nitfg, num_iter, n_points, ix, jy, nlev)
          call asad_diffun(nl)
          call asad_jac(n_points)
          call asad_impact(n_points, ix, jy, nlev)
        else
          if (method == 2) then
            cmessage(:) = 'QSSA not in UM6.5 build'
            errcode = 1
            call ereport('ASAD_CDRIVE', errcode, cmessage)
          else
            if (method == 3) then
              if (gphot) then
                call ukca_photol(prt, n_points)
              end if
              if (ukca_config%l_ukca_het_psc) then
                call ukca_hetero(n_points, have_nat, stratflag)
              end if
              call asad_ftoy(gfirst, nitfg, num_iter, n_points, ix, jy, nlev)
              call asad_spmjpdriv(ix, jy, nlev, n_points, num_iter)
              if (ukca_config%l_ukca_debug_asad) then
                if (ukca_config%l_ukca_asad_columns) then
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, "('Iterations in spmjpdriv = ',I0,' ix = ',I0,' jy = ',I0)") num_iter, ix, jy
                  call umprint(ummessage, src='asad_cdrive')
                else
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Write_Stmt
                  WRITE(umMessage, "('Iterations in spmjpdriv = ',I0,' k = ',I0)") num_iter, nlev
                  call umprint(ummessage, src='asad_cdrive')
                end if
              end if
            else
              if (method == 5) then
                if (gphot) then
                  call ukca_photol(prt, n_points)
                end if
                if (ukca_config%l_ukca_het_psc) then
                  call ukca_hetero(n_points, have_nat, stratflag)
                end if
                call asad_ftoy(gfirst, nitfg, num_iter, n_points, ix, jy, nlev)
                call asad_bedriv(ix, jy, n_points, nlev)
              else
                if (method == 10) then
                  cmessage(:) = 'NAG BDF not in UM6.5 build'
                  errcode = 1
                  call ereport('ASAD_CDRIVE', errcode, cmessage)
                else
                  if (method == 11) then
                    cmessage(:) = 'SVODE not in UM6.5 build'
                    errcode = 1
                    call ereport('ASAD_CDRIVE', errcode, cmessage)
                  end if
                end if
              end if
            end if
          end if
        end if
        if (ukca_config%l_ukca_het_psc) then
          call ukca_solidphase(n_points)
        end if
      enddo
    else
      method = 0
      if (gphot) then
        call ukca_photol(prt, n_points)
      end if
      if (ukca_config%l_ukca_het_psc) then
        call ukca_hetero(n_points, have_nat, stratflag)
      end if
      call asad_ftoy(first_call, nit0, num_iter, n_points, ix, jy, nlev)
      call asad_diffun(nl)
      if (ukca_config%l_ukca_het_psc) then
        call ukca_solidphase(n_points)
      end if
      if (ukca_config%l_ukca_het_psc) then
        call asad_posthet()
      end if
    end if
    !NO!$acc kernels
    !NO!$acc loop seq
    do js = 1, jpspec, 1
      if (ctype(js) == jpif(:)) then
        ifam = moffam(js)
        itr = madvtr(js)
        iodd = nodd(js)
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          if (linfam(jl,itr)) then
            f(jl,ifam) = f(jl,ifam) - iodd * f(jl,itr)
          end if
        enddo
      end if
    enddo
    !NO!$acc loop seq
    do jtr = 1, jpcspf, 1
      if (method /= 0) then
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          cdot(jl,jtr) = (f(jl,jtr) - ftr(jl,jtr)) / (cdt * ncsteps)
          ftr(jl,jtr) = f(jl,jtr)
        enddo
      else
        !NO!$acc loop gang vector independent
        do jl = 1, n_points, 1
          cdot(jl,jtr) = fdot(jl,jtr)
          ftr(jl,jtr) = f(jl,jtr)
        enddo
      end if
    enddo
    !NO!$acc end kernels
    if (lvmr) then
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do jtr = 1, jpcspf, 1
        do jl = 1, n_points, 1
          ftr(jl,jtr) = ftr(jl,jtr) / tnd(jl)
          cdot(jl,jtr) = cdot(jl,jtr) / tnd(jl)
        enddo
      enddo
      !NO!$acc end kernels
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_cdrive

end module asad_cdrive_mod
