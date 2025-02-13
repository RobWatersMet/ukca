module ukca_conserve_mod
  use ukca_config_specification_mod, only : ukca_config
  use yomhook, only : dr_hook, lhook
  use nvtx
  use parkind1, only : jpim, jprb
  use umPrintMgr, only : printstatus, prstatus_normal, prstatus_oper, ummessage, umprint
  implicit none
  character(len = *), parameter, private :: ModuleName = 'UKCA_CONSERVE_MOD'
  public

  contains
  subroutine ukca_conserve(row_length, rows, model_levels, ntracers, tracers, pres, drain, crain, direction)
    use asad_mod, only : advt, jpctr
    use ukca_constants, only : c_br, c_brcl, c_bro, c_brono2, c_c2h6, c_c3h8, c_ccl4, c_cf2cl2, c_cf2clbr, c_cf2clcfcl2, c_cf3br, c_cfcl3, c_ch2br2, c_ch4, c_chf2cl, c_cl, c_cl2o2, c_clo, c_clono2, c_etcho, c_etco3, c_etoo, c_etooh, c_h, c_h2, c_h2o, c_h2o2, c_hbr, c_hcho, c_hcl, c_ho2, c_ho2no2, c_hobr, c_hocl, c_hono, c_hono2, c_mebr, c_meccl3, c_mecho, c_mecl, c_mecoch2oo, c_mecoch2ooh, c_meono2, c_meoo, c_meooh, c_oclo, c_oh, c_pan, c_ppan, c_proo, c_prooh
    use ereport_mod, only : ereport
    use errormessagelength_mod, only : errormessagelength
    integer, parameter :: max_comp = 50
    real, parameter :: adjust_pres = 500.0
    real, parameter :: washout_limit = 10000.0
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: row_length
    integer, intent(in) :: rows
    integer, intent(in) :: model_levels
    integer, intent(in) :: ntracers
    logical, intent(in) :: direction
    real, dimension(row_length,rows,model_levels), intent(in) :: pres
    real, dimension(row_length,rows,model_levels), intent(in) :: drain
    real, dimension(row_length,rows,model_levels), intent(in) :: crain
    real, dimension(row_length,rows,model_levels,ntracers), intent(inout) :: tracers
    integer :: m
    integer :: i
    integer :: j
    integer :: k
    integer, save :: n_toth = 0
    integer, save :: n_n2o = 0
    logical, save :: firstcall = .true.
    integer, save :: ncl_tracers
    integer, save :: nbr_tracers
    integer, save :: nh_tracers
    real, allocatable, dimension(:,:,:), save :: total_cl
    real, allocatable, dimension(:,:,:), save :: total_br
    real, allocatable, dimension(:,:,:), save :: total_h
    integer, dimension(max_comp), save :: cl_tracers
    integer, dimension(max_comp), save :: br_tracers
    integer, dimension(max_comp), save :: h_tracers
    real, dimension(max_comp), save :: c_cl_tracers
    real, dimension(max_comp), save :: c_br_tracers
    real, dimension(max_comp), save :: c_h_tracers
    real, dimension(max_comp), save :: cl_validity
    real, dimension(max_comp), save :: br_validity
    real, dimension(max_comp), save :: h_validity
    logical, dimension(max_comp), save :: contains_bromine
    logical, dimension(max_comp), save :: do_not_change
    real, allocatable, dimension(:,:,:) :: corrfac
    integer :: errcode
    CHARACTER(LEN = errormessagelength) :: cmessage
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'UKCA_CONSERVE'
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
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19
    integer :: idx_20
    integer :: idx_21
    integer :: idx_22
    integer :: idx_23
    integer :: idx_24
    integer :: idx_25
    integer :: idx_26
    integer :: idx_27
    integer :: idx_28
    integer :: idx_29
    integer :: idx_30
    integer :: idx_31
    integer :: idx_32
    integer :: idx_33
    integer :: idx_34
    integer :: idx_35
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39
    integer :: idx_40
    integer :: idx_41
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50
    integer :: idx_51
    integer :: idx_52
    integer :: idx_53
    integer :: idx_54
    integer :: idx_55
    integer :: idx_56
    integer :: idx_57
    integer :: idx_58
    integer :: idx_59
    integer :: idx_60
    integer :: idx_61
    integer :: idx_62
    integer :: idx_63
    integer :: idx_64
    integer :: idx_65
    integer :: idx_66
    integer :: widx3
    integer :: widx2
    integer :: widx1
    integer :: widx3_1
    integer :: widx2_1
    integer :: widx1_1
    integer :: widx3_2
    integer :: widx2_2
    integer :: widx1_2

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (firstcall) then
      if (ntracers /= jpctr) then
        errcode = 1
        cmessage(:) = ' Number of tracers does not agree with jpctr'
        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported statement: Write_Stmt
        WRITE(umMessage, '(1X,A,1X,I4,1X,I4)') cmessage, ntracers, jpctr
        call umprint(ummessage, src='ukca_conserve_mod')
        call ereport('UKCA_CONSERVE_MOD:UKCA_CONSERVE', errcode, cmessage)
      end if
      !$acc kernels
      ncl_tracers = 0
      nbr_tracers = 0
      nh_tracers = 0
      !$acc loop gang vector independent
      do idx = 1, max_comp, 1
        contains_bromine(idx) = .false.
      enddo
      !$acc loop gang vector independent
      do idx_1 = 1, max_comp, 1
        do_not_change(idx_1) = .false.
      enddo
      !$acc loop gang vector independent
      do idx_2 = 1, max_comp, 1
        cl_tracers(idx_2) = 0
      enddo
      !$acc loop gang vector independent
      do idx_3 = 1, max_comp, 1
        br_tracers(idx_3) = 0
      enddo
      !$acc loop gang vector independent
      do idx_4 = 1, max_comp, 1
        h_tracers(idx_4) = 0
      enddo
      !$acc loop gang vector independent
      do idx_5 = 1, max_comp, 1
        cl_validity(idx_5) = 1.0
      enddo
      !$acc loop gang vector independent
      do idx_6 = 1, max_comp, 1
        br_validity(idx_6) = 1.0
      enddo
      !$acc loop gang vector independent
      do idx_7 = 1, max_comp, 1
        h_validity(idx_7) = 1.0
      enddo
      !$acc loop gang vector independent
      do idx_8 = 1, max_comp, 1
        c_br_tracers(idx_8) = 0.0
      enddo
      !$acc loop gang vector independent
      do idx_9 = 1, max_comp, 1
        c_cl_tracers(idx_9) = 0.0
      enddo
      !$acc loop gang vector independent
      do idx_10 = 1, max_comp, 1
        c_h_tracers(idx_10) = 0.0
      enddo
      !$acc end kernels
      do m = 1, jpctr, 1
        if (advt(m) == 'CF2Cl2    ') then
          ncl_tracers = ncl_tracers + 1
          cl_tracers(ncl_tracers) = m
          c_cl_tracers(ncl_tracers) = c_cf2cl2
          cl_validity(ncl_tracers) = 2.0
        else
          if (advt(m) == 'CFCl3     ') then
            ncl_tracers = ncl_tracers + 1
            cl_tracers(ncl_tracers) = m
            c_cl_tracers(ncl_tracers) = c_cfcl3
            cl_validity(ncl_tracers) = 3.0
          else
            if (advt(m) == 'Clx       ' .OR. advt(m) == 'ClO       ') then
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_clo
            else
              if (advt(m) == 'Cl        ') then
                ncl_tracers = ncl_tracers + 1
                cl_tracers(ncl_tracers) = m
                c_cl_tracers(ncl_tracers) = c_cl
              else
                if (advt(m) == 'Cl2O2     ') then
                  ncl_tracers = ncl_tracers + 1
                  cl_tracers(ncl_tracers) = m
                  c_cl_tracers(ncl_tracers) = c_cl2o2
                  cl_validity(ncl_tracers) = 2.0
                else
                  if (advt(m) == 'HCl       ') then
                    ncl_tracers = ncl_tracers + 1
                    cl_tracers(ncl_tracers) = m
                    c_cl_tracers(ncl_tracers) = c_hcl
                    nh_tracers = nh_tracers + 1
                    h_tracers(nh_tracers) = m
                    c_h_tracers(nh_tracers) = c_hcl
                    do_not_change(nh_tracers) = .true.
                  else
                    if (advt(m) == 'HOCl      ') then
                      ncl_tracers = ncl_tracers + 1
                      cl_tracers(ncl_tracers) = m
                      c_cl_tracers(ncl_tracers) = c_hocl
                      nh_tracers = nh_tracers + 1
                      h_tracers(nh_tracers) = m
                      c_h_tracers(nh_tracers) = c_hocl
                      do_not_change(nh_tracers) = .true.
                    else
                      if (advt(m) == 'OClO      ') then
                        ncl_tracers = ncl_tracers + 1
                        cl_tracers(ncl_tracers) = m
                        c_cl_tracers(ncl_tracers) = c_oclo
                      else
                        if (advt(m) == 'BrCl      ') then
                          ncl_tracers = ncl_tracers + 1
                          cl_tracers(ncl_tracers) = m
                          c_cl_tracers(ncl_tracers) = c_brcl
                          contains_bromine(ncl_tracers) = .true.
                          nbr_tracers = nbr_tracers + 1
                          br_tracers(nbr_tracers) = m
                          c_br_tracers(nbr_tracers) = c_brcl
                        else
                          if (advt(m) == 'ClONO2    ') then
                            ncl_tracers = ncl_tracers + 1
                            cl_tracers(ncl_tracers) = m
                            c_cl_tracers(ncl_tracers) = c_clono2
                          else
                            if (advt(m) == 'CF2ClCFCl2') then
                              ncl_tracers = ncl_tracers + 1
                              cl_tracers(ncl_tracers) = m
                              c_cl_tracers(ncl_tracers) = c_cf2clcfcl2
                              cl_validity(ncl_tracers) = 3.0
                            else
                              if (advt(m) == 'CHF2Cl    ') then
                                ncl_tracers = ncl_tracers + 1
                                cl_tracers(ncl_tracers) = m
                                c_cl_tracers(ncl_tracers) = c_chf2cl
                                nh_tracers = nh_tracers + 1
                                h_tracers(nh_tracers) = m
                                c_h_tracers(nh_tracers) = c_chf2cl
                                do_not_change(nh_tracers) = .true.
                              else
                                if (advt(m) == 'MeCCl3    ') then
                                  ncl_tracers = ncl_tracers + 1
                                  cl_tracers(ncl_tracers) = m
                                  c_cl_tracers(ncl_tracers) = c_meccl3
                                  cl_validity(ncl_tracers) = 3.0
                                  nh_tracers = nh_tracers + 1
                                  h_tracers(nh_tracers) = m
                                  c_h_tracers(nh_tracers) = c_meccl3
                                  do_not_change(nh_tracers) = .true.
                                else
                                  if (advt(m) == 'CCl4      ') then
                                    ncl_tracers = ncl_tracers + 1
                                    cl_tracers(ncl_tracers) = m
                                    c_cl_tracers(ncl_tracers) = c_ccl4
                                    cl_validity(ncl_tracers) = 4.0
                                  else
                                    if (advt(m) == 'MeCl      ') then
                                      ncl_tracers = ncl_tracers + 1
                                      cl_tracers(ncl_tracers) = m
                                      c_cl_tracers(ncl_tracers) = c_mecl
                                      nh_tracers = nh_tracers + 1
                                      h_tracers(nh_tracers) = m
                                      c_h_tracers(nh_tracers) = c_mecl
                                      do_not_change(nh_tracers) = .true.
                                    else
                                      if (advt(m) == 'CF2ClBr   ') then
                                        ncl_tracers = ncl_tracers + 1
                                        cl_tracers(ncl_tracers) = m
                                        c_cl_tracers(ncl_tracers) = c_cf2clbr
                                        contains_bromine(ncl_tracers) = .true.
                                        nbr_tracers = nbr_tracers + 1
                                        br_tracers(nbr_tracers) = m
                                        c_br_tracers(nbr_tracers) = c_cf2clbr
                                      else
                                        if (advt(m) == 'Brx       ' .OR. advt(m) == 'BrO       ') then
                                          nbr_tracers = nbr_tracers + 1
                                          br_tracers(nbr_tracers) = m
                                          c_br_tracers(nbr_tracers) = c_bro
                                        else
                                          if (advt(m) == 'Br        ') then
                                            nbr_tracers = nbr_tracers + 1
                                            br_tracers(nbr_tracers) = m
                                            c_br_tracers(nbr_tracers) = c_br
                                          else
                                            if (advt(m) == 'HOBr      ') then
                                              nbr_tracers = nbr_tracers + 1
                                              br_tracers(nbr_tracers) = m
                                              c_br_tracers(nbr_tracers) = c_hobr
                                              nh_tracers = nh_tracers + 1
                                              h_tracers(nh_tracers) = m
                                              c_h_tracers(nh_tracers) = c_hobr
                                              do_not_change(nh_tracers) = .true.
                                            else
                                              if (advt(m) == 'BrONO2    ') then
                                                nbr_tracers = nbr_tracers + 1
                                                br_tracers(nbr_tracers) = m
                                                c_br_tracers(nbr_tracers) = c_brono2
                                              else
                                                if (advt(m) == 'MeBr      ') then
                                                  nbr_tracers = nbr_tracers + 1
                                                  br_tracers(nbr_tracers) = m
                                                  c_br_tracers(nbr_tracers) = c_mebr
                                                  nh_tracers = nh_tracers + 1
                                                  h_tracers(nh_tracers) = m
                                                  c_h_tracers(nh_tracers) = c_mebr
                                                  do_not_change(nh_tracers) = .true.
                                                else
                                                  if (advt(m) == 'HBr       ') then
                                                    nbr_tracers = nbr_tracers + 1
                                                    br_tracers(nbr_tracers) = m
                                                    c_br_tracers(nbr_tracers) = c_hbr
                                                    nh_tracers = nh_tracers + 1
                                                    h_tracers(nh_tracers) = m
                                                    c_h_tracers(nh_tracers) = c_hbr
                                                    do_not_change(nh_tracers) = .true.
                                                  else
                                                    if (advt(m) == 'CF3Br     ') then
                                                      nbr_tracers = nbr_tracers + 1
                                                      br_tracers(nbr_tracers) = m
                                                      c_br_tracers(nbr_tracers) = c_cf3br
                                                    else
                                                      if (advt(m) == 'CH2Br2    ') then
                                                        nbr_tracers = nbr_tracers + 1
                                                        br_tracers(nbr_tracers) = m
                                                        c_br_tracers(nbr_tracers) = c_ch2br2
                                                        br_validity(nbr_tracers) = 2.0
                                                        nh_tracers = nh_tracers + 1
                                                        h_tracers(nh_tracers) = m
                                                        c_h_tracers(nh_tracers) = c_ch2br2
                                                        h_validity(nh_tracers) = 2.0
                                                        do_not_change(nh_tracers) = .true.
                                                      else
                                                        if (advt(m) == 'TOTH      ') then
                                                          n_toth = m
                                                        else
                                                          if (advt(m) == 'H2O       ') then
                                                            nh_tracers = nh_tracers + 1
                                                            h_tracers(nh_tracers) = m
                                                            c_h_tracers(nh_tracers) = c_h2o
                                                            h_validity(nh_tracers) = 2.0
                                                          else
                                                            if (advt(m) == 'CH4       ') then
                                                              nh_tracers = nh_tracers + 1
                                                              h_tracers(nh_tracers) = m
                                                              c_h_tracers(nh_tracers) = c_ch4
                                                              h_validity(nh_tracers) = 4.0
                                                              do_not_change(nh_tracers) = .true.
                                                            else
                                                              if (advt(m) == 'H2        ') then
                                                                nh_tracers = nh_tracers + 1
                                                                h_tracers(nh_tracers) = m
                                                                c_h_tracers(nh_tracers) = c_h2
                                                                h_validity(nh_tracers) = 2.0
                                                              else
                                                                if (advt(m) == 'HO2NO2    ') then
                                                                  nh_tracers = nh_tracers + 1
                                                                  h_tracers(nh_tracers) = m
                                                                  c_h_tracers(nh_tracers) = c_ho2no2
                                                                  do_not_change(nh_tracers) = .true.
                                                                else
                                                                  if (advt(m) == 'HONO2     ') then
                                                                    nh_tracers = nh_tracers + 1
                                                                    h_tracers(nh_tracers) = m
                                                                    c_h_tracers(nh_tracers) = c_hono2
                                                                    do_not_change(nh_tracers) = .true.
                                                                  else
                                                                    if (advt(m) == 'H2O2      ') then
                                                                      nh_tracers = nh_tracers + 1
                                                                      h_tracers(nh_tracers) = m
                                                                      c_h_tracers(nh_tracers) = c_h2o2
                                                                      h_validity(nh_tracers) = 2.0
                                                                    else
                                                                      if (advt(m) == 'HCHO      ') then
                                                                        nh_tracers = nh_tracers + 1
                                                                        h_tracers(nh_tracers) = m
                                                                        c_h_tracers(nh_tracers) = c_hcho
                                                                        h_validity(nh_tracers) = 2.0
                                                                      else
                                                                        if (advt(m) == 'MeOOH     ') then
                                                                          nh_tracers = nh_tracers + 1
                                                                          h_tracers(nh_tracers) = m
                                                                          c_h_tracers(nh_tracers) = c_meooh
                                                                          h_validity(nh_tracers) = 4.0
                                                                        else
                                                                          if (advt(m) == 'HONO      ') then
                                                                            nh_tracers = nh_tracers + 1
                                                                            h_tracers(nh_tracers) = m
                                                                            c_h_tracers(nh_tracers) = c_hono
                                                                            do_not_change(nh_tracers) = .true.
                                                                          else
                                                                            if (advt(m) == 'C2H6      ') then
                                                                              nh_tracers = nh_tracers + 1
                                                                              h_tracers(nh_tracers) = m
                                                                              c_h_tracers(nh_tracers) = c_c2h6
                                                                              h_validity(nh_tracers) = 6.0
                                                                            else
                                                                              if (advt(m) == 'EtOOH     ') then
                                                                                nh_tracers = nh_tracers + 1
                                                                                h_tracers(nh_tracers) = m
                                                                                c_h_tracers(nh_tracers) = c_etooh
                                                                                h_validity(nh_tracers) = 6.0
                                                                              else
                                                                                if (advt(m) == 'MeCHO     ') then
                                                                                  nh_tracers = nh_tracers + 1
                                                                                  h_tracers(nh_tracers) = m
                                                                                  c_h_tracers(nh_tracers) = c_mecho
                                                                                  h_validity(nh_tracers) = 4.0
                                                                                else
                                                                                  if (advt(m) == 'PAN       ') then
                                                                                    nh_tracers = nh_tracers + 1
                                                                                    h_tracers(nh_tracers) = m
                                                                                    c_h_tracers(nh_tracers) = c_pan
                                                                                    h_validity(nh_tracers) = 3.0
                                                                                    do_not_change(nh_tracers) = .true.
                                                                                  else
                                                                                    if (advt(m) == 'C3H8      ') then
                                                                                      nh_tracers = nh_tracers + 1
                                                                                      h_tracers(nh_tracers) = m
                                                                                      c_h_tracers(nh_tracers) = c_c3h8
                                                                                      h_validity(nh_tracers) = 8.0
                                                                                    else
                                                                                      if (advt(m) == 'n-PrOOH   ' .OR. advt(m) == 'i-PrOOH   ') then
                                                                                        nh_tracers = nh_tracers + 1
                                                                                        h_tracers(nh_tracers) = m
                                                                                        c_h_tracers(nh_tracers) = c_prooh
                                                                                        h_validity(nh_tracers) = 8.0
                                                                                      else
                                                                                        if (advt(m) == 'EtCHO     ' .OR. advt(m) == 'Me2CO     ') then
                                                                                          nh_tracers = nh_tracers + 1
                                                                                          h_tracers(nh_tracers) = m
                                                                                          c_h_tracers(nh_tracers) = c_etcho
                                                                                          h_validity(nh_tracers) = 6.0
                                                                                        else
                                                                                          if (advt(m) == 'MeCOCH2OOH') then
                                                                                            nh_tracers = nh_tracers + 1
                                                                                            h_tracers(nh_tracers) = m
                                                                                            c_h_tracers(nh_tracers) = c_mecoch2ooh
                                                                                            h_validity(nh_tracers) = 6.0
                                                                                          else
                                                                                            if (advt(m) == 'PPAN      ') then
                                                                                              nh_tracers = nh_tracers + 1
                                                                                              h_tracers(nh_tracers) = m
                                                                                              c_h_tracers(nh_tracers) = c_ppan
                                                                                              h_validity(nh_tracers) = 5.0
                                                                                              do_not_change(nh_tracers) = .true.
                                                                                            else
                                                                                              if (advt(m) == 'MeONO2    ') then
                                                                                                nh_tracers = nh_tracers + 1
                                                                                                h_tracers(nh_tracers) = m
                                                                                                c_h_tracers(nh_tracers) = c_meono2
                                                                                                h_validity(nh_tracers) = 3.0
                                                                                                do_not_change(nh_tracers) = .true.
                                                                                              else
                                                                                                if (advt(m) == 'H         ') then
                                                                                                  nh_tracers = nh_tracers + 1
                                                                                                  h_tracers(nh_tracers) = m
                                                                                                  c_h_tracers(nh_tracers) = c_h
                                                                                                else
                                                                                                  if (advt(m) == 'OH        ') then
                                                                                                    nh_tracers = nh_tracers + 1
                                                                                                    h_tracers(nh_tracers) = m
                                                                                                    c_h_tracers(nh_tracers) = c_oh
                                                                                                  else
                                                                                                    if (advt(m) == 'HO2       ' .OR. advt(m) == 'HOx       ') then
                                                                                                      nh_tracers = nh_tracers + 1
                                                                                                      h_tracers(nh_tracers) = m
                                                                                                      c_h_tracers(nh_tracers) = c_ho2
                                                                                                    else
                                                                                                      if (advt(m) == 'MeOO      ') then
                                                                                                        nh_tracers = nh_tracers + 1
                                                                                                        h_tracers(nh_tracers) = m
                                                                                                        c_h_tracers(nh_tracers) = c_meoo
                                                                                                        h_validity(nh_tracers) = 3.0
                                                                                                      else
                                                                                                        if (advt(m) == 'EtOO      ') then
                                                                                                          nh_tracers = nh_tracers + 1
                                                                                                          h_tracers(nh_tracers) = m
                                                                                                          c_h_tracers(nh_tracers) = c_etoo
                                                                                                          h_validity(nh_tracers) = 5.0
                                                                                                        else
                                                                                                          if (advt(m) == 'i-PrOO    ' .OR. advt(m) == 'n-PrOO    ') then
                                                                                                            nh_tracers = nh_tracers + 1
                                                                                                            h_tracers(nh_tracers) = m
                                                                                                            c_h_tracers(nh_tracers) = c_proo
                                                                                                            h_validity(nh_tracers) = 7.0
                                                                                                          else
                                                                                                            if (advt(m) == 'EtCO3     ') then
                                                                                                              nh_tracers = nh_tracers + 1
                                                                                                              h_tracers(nh_tracers) = m
                                                                                                              c_h_tracers(nh_tracers) = c_etco3
                                                                                                              h_validity(nh_tracers) = 5.0
                                                                                                            else
                                                                                                              if (advt(m) == 'MeCOCH2OO ') then
                                                                                                                nh_tracers = nh_tracers + 1
                                                                                                                h_tracers(nh_tracers) = m
                                                                                                                c_h_tracers(nh_tracers) = c_mecoch2oo
                                                                                                                h_validity(nh_tracers) = 5.0
                                                                                                              else
                                                                                                                if (advt(m) == 'Ox        ' .OR. advt(m) == 'O3        ') then
                                                                                                                else
                                                                                                                  if (advt(m) == 'N2O       ') then
                                                                                                                    n_n2o = m
                                                                                                                  else
                                                                                                                    if (printstatus >= prstatus_oper) then
                                                                                                                      cmessage(:) = ' CONSERVE: SPECIES ' // advt(m) // ' not treated in CASE'
                                                                                                                      errcode = -1 * m
                                                                                                                      call ereport('UKCA_CONSERVE_MOD:UKCA_CONSERVE', errcode, cmessage)
                                                                                                                    end if
                                                                                                                  end if
                                                                                                                end if
                                                                                                              end if
                                                                                                            end if
                                                                                                          end if
                                                                                                        end if
                                                                                                      end if
                                                                                                    end if
                                                                                                  end if
                                                                                                end if
                                                                                              end if
                                                                                            end if
                                                                                          end if
                                                                                        end if
                                                                                      end if
                                                                                    end if
                                                                                  end if
                                                                                end if
                                                                              end if
                                                                            end if
                                                                          end if
                                                                        end if
                                                                      end if
                                                                    end if
                                                                  end if
                                                                end if
                                                              end if
                                                            end if
                                                          end if
                                                        end if
                                                      end if
                                                    end if
                                                  end if
                                                end if
                                              end if
                                            end if
                                          end if
                                        end if
                                      end if
                                    end if
                                  end if
                                end if
                              end if
                            end if
                          end if
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
      firstcall = .false.
    end if
    if (direction) then
      if (.NOT.ALLOCATED(total_br)) then
        ALLOCATE(total_br(1:row_length,1:rows,1:model_levels))
      end if
      if (.NOT.ALLOCATED(total_cl)) then
        ALLOCATE(total_cl(1:row_length,1:rows,1:model_levels))
      end if
      !$acc kernels
      !$acc loop gang vector independent collapse(3)
      do idx_11 = LBOUND(total_br, dim=3), UBOUND(total_br, dim=3), 1
        do idx_12 = LBOUND(total_br, dim=2), UBOUND(total_br, dim=2), 1
          do idx_13 = LBOUND(total_br, dim=1), UBOUND(total_br, dim=1), 1
            total_br(idx_13,idx_12,idx_11) = 0.0
          enddo
        enddo
      enddo
      !$acc loop gang vector independent collapse(3)
      do idx_14 = LBOUND(total_cl, dim=3), UBOUND(total_cl, dim=3), 1
        do idx_15 = LBOUND(total_cl, dim=2), UBOUND(total_cl, dim=2), 1
          do idx_16 = LBOUND(total_cl, dim=1), UBOUND(total_cl, dim=1), 1
            total_cl(idx_16,idx_15,idx_14) = 0.0
          enddo
        enddo
      enddo
      !$acc loop seq
      do m = 1, nbr_tracers, 1
        !$acc loop gang vector independent collapse(3)
        do idx_17 = LBOUND(total_br, dim=3), UBOUND(total_br, dim=3), 1
          do idx_18 = LBOUND(total_br, dim=2), UBOUND(total_br, dim=2), 1
            do idx_19 = LBOUND(total_br, dim=1), UBOUND(total_br, dim=1), 1
              total_br(idx_19,idx_18,idx_17) = total_br(idx_19,idx_18,idx_17) + tracers(idx_19 + (LBOUND(tracers, dim=1) - LBOUND(total_br, dim=1)),idx_18 + (LBOUND(tracers, dim=2) - LBOUND(total_br, dim=2)),idx_17 + (LBOUND(tracers, dim=3) - LBOUND(total_br, dim=3)),br_tracers(m)) * br_validity(m) / c_br_tracers(m)
            enddo
          enddo
        enddo
      enddo
      !$acc loop seq
      do m = 1, ncl_tracers, 1
        !$acc loop gang vector independent collapse(3)
        do idx_20 = LBOUND(total_cl, dim=3), UBOUND(total_cl, dim=3), 1
          do idx_21 = LBOUND(total_cl, dim=2), UBOUND(total_cl, dim=2), 1
            do idx_22 = LBOUND(total_cl, dim=1), UBOUND(total_cl, dim=1), 1
              total_cl(idx_22,idx_21,idx_20) = total_cl(idx_22,idx_21,idx_20) + tracers(idx_22 + (LBOUND(tracers, dim=1) - LBOUND(total_cl, dim=1)),idx_21 + (LBOUND(tracers, dim=2) - LBOUND(total_cl, dim=2)),idx_20 + (LBOUND(tracers, dim=3) - LBOUND(total_cl, dim=3)),cl_tracers(m)) * cl_validity(m) / c_cl_tracers(m)
            enddo
          enddo
        enddo
      enddo
      !$acc end kernels
      if (ukca_config%l_ukca_h2o_feedback .AND. ukca_config%l_ukca_conserve_h) then
        ALLOCATE(total_h(1:row_length,1:rows,1:model_levels))
        !$acc kernels
        !$acc loop gang vector independent collapse(3)
        do idx_23 = LBOUND(total_h, dim=3), UBOUND(total_h, dim=3), 1
          do idx_24 = LBOUND(total_h, dim=2), UBOUND(total_h, dim=2), 1
            do idx_25 = LBOUND(total_h, dim=1), UBOUND(total_h, dim=1), 1
              total_h(idx_25,idx_24,idx_23) = 0.0
            enddo
          enddo
        enddo
        !$acc end kernels
        !$acc kernels
        !$acc loop seq
        do m = 1, nh_tracers, 1
          !$acc loop gang vector independent collapse(3)
          do idx_26 = LBOUND(total_h, dim=3), UBOUND(total_h, dim=3), 1
            do idx_27 = LBOUND(total_h, dim=2), UBOUND(total_h, dim=2), 1
              do idx_28 = LBOUND(total_h, dim=1), UBOUND(total_h, dim=1), 1
                total_h(idx_28,idx_27,idx_26) = total_h(idx_28,idx_27,idx_26) + tracers(idx_28 + (LBOUND(tracers, dim=1) - LBOUND(total_h, dim=1)),idx_27 + (LBOUND(tracers, dim=2) - LBOUND(total_h, dim=2)),idx_26 + (LBOUND(tracers, dim=3) - LBOUND(total_h, dim=3)),h_tracers(m)) * h_validity(m) / c_h_tracers(m)
              enddo
            enddo
          enddo
        enddo
        !$acc end kernels
        if (n_toth > 0) then
          !$acc kernels
          !$acc loop gang vector independent collapse(3)
          do widx3 = 1, model_levels, 1
            do widx2 = 1, rows, 1
              do widx1 = 1, row_length, 1
                if (pres(widx1,widx2,widx3) > adjust_pres) then
                  tracers(widx1,widx2,widx3,n_toth) = total_h(LBOUND(total_h, dim=1) + widx1 - 1,LBOUND(total_h, dim=2) + widx2 - 1,LBOUND(total_h, dim=3) + widx3 - 1)
                end if
              enddo
            enddo
          enddo
          !$acc loop gang vector independent collapse(3)
          do widx3_1 = 1, model_levels, 1
            do widx2_1 = 1, rows, 1
              do widx1_1 = 1, row_length, 1
                if (pres(widx1_1,widx2_1,widx3_1) <= adjust_pres) then
                  total_h(LBOUND(total_h, dim=1) + widx1_1 - 1,LBOUND(total_h, dim=2) + widx2_1 - 1,LBOUND(total_h, dim=3) + widx3_1 - 1) = tracers(widx1_1,widx2_1,widx3_1,n_toth)
                end if
              enddo
            enddo
          enddo
          !$acc end kernels
        end if
      end if
    else
      if (n_n2o > 0 .AND. model_levels /= 1) then
        !$acc kernels
        !$acc loop gang vector independent collapse(2)
        do idx_29 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
          do idx_30 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
            tracers(idx_30,idx_29,model_levels,n_n2o) = 0.0
          enddo
        enddo
        !$acc end kernels
      end if
      ALLOCATE(corrfac(1:row_length,1:rows,1:model_levels))
      if (nbr_tracers > 0) then
        !$acc kernels
        !$acc loop gang vector independent collapse(3)
        do idx_31 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_32 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_33 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_33,idx_32,idx_31) = 0.0
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, nbr_tracers, 1
          !$acc loop gang vector independent collapse(3)
          do idx_34 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
            do idx_35 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
              do idx_36 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
                corrfac(idx_36,idx_35,idx_34) = corrfac(idx_36,idx_35,idx_34) + tracers(idx_36 + (LBOUND(tracers, dim=1) - LBOUND(corrfac, dim=1)),idx_35 + (LBOUND(tracers, dim=2) - LBOUND(corrfac, dim=2)),idx_34 + (LBOUND(tracers, dim=3) - LBOUND(corrfac, dim=3)),br_tracers(m)) * br_validity(m) / c_br_tracers(m)
              enddo
            enddo
          enddo
        enddo
        !$acc loop gang vector independent collapse(3)
        do idx_37 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_38 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_39 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_39,idx_38,idx_37) = total_br(idx_39 + (LBOUND(total_br, dim=1) - LBOUND(corrfac, dim=1)),idx_38 + (LBOUND(total_br, dim=2) - LBOUND(corrfac, dim=2)),idx_37 + (LBOUND(total_br, dim=3) - LBOUND(corrfac, dim=3))) / corrfac(idx_39,idx_38,idx_37)
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, nbr_tracers, 1
          !$acc loop gang vector independent collapse(3)
          do i = 1, row_length, 1
            do j = 1, rows, 1
              do k = 2, model_levels, 1
                if (drain(i,j,k) + crain(i,j,k) == 0.0) then
                  tracers(i,j,k,br_tracers(m)) = tracers(i,j,k,br_tracers(m)) * corrfac(i,j,k)
                else
                  corrfac(i,j,k) = 1.0
                end if
              enddo
            enddo
          enddo
        enddo
        !$acc end kernels
        if ((MINVAL(corrfac) < 0.9 .OR. MAXVAL(corrfac) > 1.1) .AND. printstatus >= prstatus_normal) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,2F8.4)') 'Correct bromine ', MINVAL(corrfac), MAXVAL(corrfac)
          call umprint(ummessage, src='ukca_conserve_mod')
        end if
      end if
      if (ncl_tracers > 0) then
        !$acc kernels
        !$acc loop gang vector independent collapse(3)
        do idx_40 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_41 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_42 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_42,idx_41,idx_40) = 0.0
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, ncl_tracers, 1
          if (.NOT.contains_bromine(m)) then
            !$acc loop gang vector independent collapse(3)
            do idx_43 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
              do idx_44 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
                do idx_45 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
                  corrfac(idx_45,idx_44,idx_43) = corrfac(idx_45,idx_44,idx_43) + tracers(idx_45 + (LBOUND(tracers, dim=1) - LBOUND(corrfac, dim=1)),idx_44 + (LBOUND(tracers, dim=2) - LBOUND(corrfac, dim=2)),idx_43 + (LBOUND(tracers, dim=3) - LBOUND(corrfac, dim=3)),cl_tracers(m)) * cl_validity(m) / c_cl_tracers(m)
                enddo
              enddo
            enddo
          else
            !$acc loop gang vector independent collapse(3)
            do idx_46 = LBOUND(total_cl, dim=3), UBOUND(total_cl, dim=3), 1
              do idx_47 = LBOUND(total_cl, dim=2), UBOUND(total_cl, dim=2), 1
                do idx_48 = LBOUND(total_cl, dim=1), UBOUND(total_cl, dim=1), 1
                  total_cl(idx_48,idx_47,idx_46) = total_cl(idx_48,idx_47,idx_46) - tracers(idx_48 + (LBOUND(tracers, dim=1) - LBOUND(total_cl, dim=1)),idx_47 + (LBOUND(tracers, dim=2) - LBOUND(total_cl, dim=2)),idx_46 + (LBOUND(tracers, dim=3) - LBOUND(total_cl, dim=3)),cl_tracers(m)) * cl_validity(m) / c_cl_tracers(m)
                enddo
              enddo
            enddo
          end if
        enddo
        !$acc loop gang vector independent collapse(3)
        do idx_49 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_50 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_51 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_51,idx_50,idx_49) = total_cl(idx_51 + (LBOUND(total_cl, dim=1) - LBOUND(corrfac, dim=1)),idx_50 + (LBOUND(total_cl, dim=2) - LBOUND(corrfac, dim=2)),idx_49 + (LBOUND(total_cl, dim=3) - LBOUND(corrfac, dim=3))) / corrfac(idx_51,idx_50,idx_49)
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, ncl_tracers, 1
          if (.NOT.contains_bromine(m)) then
            !$acc loop gang vector independent collapse(3)
            do i = 1, row_length, 1
              do j = 1, rows, 1
                do k = 2, model_levels, 1
                  if (crain(i,j,k) + drain(i,j,k) == 0.0) then
                    tracers(i,j,k,cl_tracers(m)) = tracers(i,j,k,cl_tracers(m)) * corrfac(i,j,k)
                  else
                    corrfac(i,j,k) = 1.0
                  end if
                enddo
              enddo
            enddo
          end if
        enddo
        !$acc end kernels
        if ((MINVAL(corrfac) < 0.9 .OR. MAXVAL(corrfac) > 1.1) .AND. printstatus >= prstatus_normal) then
          ! PSyclone CodeBlock (unsupported code) reason:
          !  - Unsupported statement: Write_Stmt
          WRITE(umMessage, '(A,2F8.4)') 'Correct chlorine ', MINVAL(corrfac), MAXVAL(corrfac)
          call umprint(ummessage, src='ukca_conserve_mod')
        end if
      end if
      if (ukca_config%l_ukca_h2o_feedback .AND. ukca_config%l_ukca_conserve_h) then
        !$acc kernels
        !$acc loop gang vector independent collapse(3)
        do idx_52 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_53 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_54 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_54,idx_53,idx_52) = 0.0
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, nh_tracers, 1
          if (.NOT.do_not_change(m)) then
            !$acc loop gang vector independent collapse(3)
            do idx_55 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
              do idx_56 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
                do idx_57 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
                  corrfac(idx_57,idx_56,idx_55) = corrfac(idx_57,idx_56,idx_55) + tracers(idx_57 + (LBOUND(tracers, dim=1) - LBOUND(corrfac, dim=1)),idx_56 + (LBOUND(tracers, dim=2) - LBOUND(corrfac, dim=2)),idx_55 + (LBOUND(tracers, dim=3) - LBOUND(corrfac, dim=3)),h_tracers(m)) * h_validity(m) / c_h_tracers(m)
                enddo
              enddo
            enddo
          else
            !$acc loop gang vector independent collapse(3)
            do idx_58 = LBOUND(total_h, dim=3), UBOUND(total_h, dim=3), 1
              do idx_59 = LBOUND(total_h, dim=2), UBOUND(total_h, dim=2), 1
                do idx_60 = LBOUND(total_h, dim=1), UBOUND(total_h, dim=1), 1
                  total_h(idx_60,idx_59,idx_58) = total_h(idx_60,idx_59,idx_58) - tracers(idx_60 + (LBOUND(tracers, dim=1) - LBOUND(total_h, dim=1)),idx_59 + (LBOUND(tracers, dim=2) - LBOUND(total_h, dim=2)),idx_58 + (LBOUND(tracers, dim=3) - LBOUND(total_h, dim=3)),h_tracers(m)) * h_validity(m) / c_h_tracers(m)
                enddo
              enddo
            enddo
          end if
        enddo
        !$acc loop gang vector independent collapse(3)
        do idx_61 = LBOUND(corrfac, dim=3), UBOUND(corrfac, dim=3), 1
          do idx_62 = LBOUND(corrfac, dim=2), UBOUND(corrfac, dim=2), 1
            do idx_63 = LBOUND(corrfac, dim=1), UBOUND(corrfac, dim=1), 1
              corrfac(idx_63,idx_62,idx_61) = total_h(idx_63 + (LBOUND(total_h, dim=1) - LBOUND(corrfac, dim=1)),idx_62 + (LBOUND(total_h, dim=2) - LBOUND(corrfac, dim=2)),idx_61 + (LBOUND(total_h, dim=3) - LBOUND(corrfac, dim=3))) / corrfac(idx_63,idx_62,idx_61)
            enddo
          enddo
        enddo
        !$acc loop gang vector independent collapse(3)
        do widx3_2 = 1, model_levels, 1
          do widx2_2 = 1, rows, 1
            do widx1_2 = 1, row_length, 1
              if (pres(widx1_2,widx2_2,widx3_2) > washout_limit) then
                corrfac(LBOUND(corrfac, dim=1) + widx1_2 - 1,LBOUND(corrfac, dim=2) + widx2_2 - 1,LBOUND(corrfac, dim=3) + widx3_2 - 1) = 1.0
              end if
            enddo
          enddo
        enddo
        !$acc loop seq
        do m = 1, nh_tracers, 1
          if (.NOT.do_not_change(m)) then
            !$acc loop gang vector independent collapse(3)
            do idx_64 = LBOUND(tracers, dim=3), UBOUND(tracers, dim=3), 1
              do idx_65 = LBOUND(tracers, dim=2), UBOUND(tracers, dim=2), 1
                do idx_66 = LBOUND(tracers, dim=1), UBOUND(tracers, dim=1), 1
                  tracers(idx_66,idx_65,idx_64,h_tracers(m)) = tracers(idx_66,idx_65,idx_64,h_tracers(m)) * corrfac(idx_66 + (LBOUND(corrfac, dim=1) - LBOUND(tracers, dim=1)),idx_65 + (LBOUND(corrfac, dim=2) - LBOUND(tracers, dim=2)),idx_64 + (LBOUND(corrfac, dim=3) - LBOUND(tracers, dim=3)))
                enddo
              enddo
            enddo
          end if
        enddo
        !$acc end kernels
      end if
      if (ALLOCATED(total_br)) then
        DEALLOCATE(total_br)
      end if
      if (ALLOCATED(total_cl)) then
        DEALLOCATE(total_cl)
      end if
      if (ALLOCATED(total_h)) then
        DEALLOCATE(total_h)
      end if
      if (ALLOCATED(corrfac)) then
        DEALLOCATE(corrfac)
      end if
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine ukca_conserve

end module ukca_conserve_mod
