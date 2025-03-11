module asad_mod
  use ukca_config_specification_mod, only : int_method_nr, ukca_config
  use yomhook, only : dr_hook, lhook
  use nvtx
  use parkind1, only : jpim, jprb
  use ukca_missing_data_mod, only : imdi
  implicit none
  real, parameter, public :: pmin = 1.0e-20
  real, parameter, public :: ptol = 1.0e-5
  real, parameter, public :: ftol = 1.0e-3
  integer, parameter, public :: kfphot = 0
  integer, parameter, public :: jpss = 16
  integer, parameter, public :: jpssr = 51
  integer, parameter, public :: nvar = 17
  integer, parameter, public :: nllv = 17
  integer, parameter, public :: ninv = 23
  integer, parameter, public :: nout = 71
  integer, parameter, public :: jpem = 9
  integer, parameter, public :: jpeq = 2
  integer, parameter, public :: jddept = 6
  integer, parameter, public :: jddepc = 5
  integer, parameter, public :: jpdwio = 56
  integer, parameter, public :: jpemio = 57
  integer, parameter, public :: jpfrpd = 100
  integer, parameter, public :: jpab = 3
  integer, parameter, public :: jpat = 7
  integer, parameter, public :: jpaj = 3
  integer, parameter, public :: jpah = 3
  integer, parameter, public :: jpspb = 6
  integer, parameter, public :: jpspt = 4
  integer, parameter, public :: jpspj = 6
  integer, parameter, public :: jpsph = 6
  integer, parameter, public :: jpmsp = jpspb
  integer, parameter, public :: jppjac = 10
  integer, parameter, public :: jpkargs = 10
  integer, parameter, public :: jprargs = 10
  integer, parameter, public :: jpcargs = 1
  integer, parameter, public :: jpcio = 55
  logical, parameter, public :: lvmr = .true.
  integer, parameter, public :: nss_o1d = 1
  integer, parameter, public :: nss_o3p = 2
  integer, parameter, public :: nss_n = 3
  integer, parameter, public :: nss_h = 4
  real, parameter, public :: tslimit = 1200.0
  integer, parameter, public :: nrsteps_max = 200
  integer, save, public :: jpctr = 0
  integer, save, public :: jpspec = 0
  integer, save, public :: jpro2 = 0
  integer, save, public :: jpcspf = 0
  integer, save, public :: jpbk = 0
  integer, save, public :: jptk = 0
  integer, save, public :: jppj = 0
  integer, save, public :: jphk = 0
  integer, save, public :: jpnr = 0
  integer, save, public :: jpdd = 0
  integer, save, public :: jpdw = 0
  real, allocatable, dimension(:), public :: wp
  real, allocatable, dimension(:), public :: co2
  real, allocatable, dimension(:,:), public :: dpd
  real, allocatable, dimension(:,:), public :: dpw
  real, allocatable, dimension(:,:), public :: emr
  real, allocatable, dimension(:,:,:), public :: fj
  real, allocatable, dimension(:,:), public :: qa
  real, allocatable, dimension(:,:), public :: ratio
  real, allocatable, dimension(:), public :: p
  real, allocatable, dimension(:), public :: t
  real, allocatable, dimension(:), public :: t300
  real, allocatable, dimension(:), public :: tnd
  real, allocatable, dimension(:), public :: pmintnd
  real, allocatable, dimension(:,:), public :: f
  real, allocatable, dimension(:,:), public :: fdot
  real, target, allocatable, public :: pd(:, :)
  real, pointer, public :: prod(:, :)
  real, pointer, public :: slos(:, :)
  real, allocatable, dimension(:,:), public :: y
  real, allocatable, dimension(:,:), public :: ydot
  real, allocatable, dimension(:,:), public :: ftilde
  real, allocatable, dimension(:,:), public :: ej
  real, allocatable, dimension(:,:), public :: rk
  real, allocatable, dimension(:,:), public :: prk
  real, allocatable, dimension(:,:,:), public :: deriv
  real, allocatable, dimension(:), public :: za
  real, allocatable, dimension(:), public :: co3
  real, allocatable, dimension(:), public :: lati
  real, allocatable, dimension(:), public :: sphno3
  real, allocatable, dimension(:), public :: sph2o
  real, allocatable, dimension(:,:,:), public :: depvel
  real, allocatable, dimension(:), public :: k298
  real, allocatable, dimension(:), public :: dhr
  real, allocatable, dimension(:,:), public :: kd298
  real, allocatable, dimension(:,:), public :: ddhr
  real, allocatable, dimension(:), public :: ct_k298
  real, allocatable, dimension(:), public :: ct_dhr
  real, allocatable, dimension(:,:), public :: ct_kd298
  real, allocatable, dimension(:,:), public :: ct_ddhr
  real, allocatable, dimension(:,:), public :: ab
  real, allocatable, dimension(:,:), public :: at
  real, allocatable, dimension(:,:), public :: aj
  real, allocatable, dimension(:,:), public :: ah
  real, allocatable, dimension(:,:), public :: ztabpd
  real, allocatable, dimension(:), public :: shno3
  real, allocatable, dimension(:), public :: sh2o
  real, allocatable, dimension(:), public :: fpsc1
  real, allocatable, dimension(:), public :: fpsc2
  real, allocatable, dimension(:,:), public :: spfj
  integer, allocatable, dimension(:), public :: madvtr
  integer, allocatable, dimension(:), public :: majors
  integer, allocatable, dimension(:), public :: moffam
  integer, allocatable, dimension(:), public :: nodd
  integer, allocatable, dimension(:), public :: nltrf
  integer, allocatable, dimension(:), public :: nltr3
  integer, allocatable, dimension(:), public :: nlnaro2
  integer, allocatable, dimension(:), public :: nlfro2
  integer, allocatable, dimension(:,:), public :: ipa
  integer, allocatable, dimension(:), public :: ipa2
  integer, allocatable, dimension(:,:), public :: nltrim
  integer, allocatable, dimension(:,:), public :: nlpdv
  integer, allocatable, dimension(:), public :: nfrpx
  integer, allocatable, dimension(:,:), public :: ntabfp
  integer, allocatable, dimension(:,:), public :: ntabpd
  integer, allocatable, dimension(:,:), public :: npdfr
  integer, allocatable, dimension(:,:), public :: ngrp
  integer, allocatable, dimension(:,:), public :: njcgrp
  integer, allocatable, dimension(:,:,:), public :: nprdx3
  integer, allocatable, dimension(:,:), public :: nprdx2
  integer, allocatable, dimension(:), public :: nprdx1
  integer, allocatable, dimension(:,:,:), public :: njacx3
  integer, allocatable, dimension(:,:), public :: njacx2
  integer, allocatable, dimension(:), public :: njacx1
  integer, allocatable, dimension(:), public :: nmpjac
  integer, allocatable, dimension(:,:), public :: npjac1
  integer, allocatable, dimension(:), public :: nbrkx
  integer, allocatable, dimension(:), public :: ntrkx
  integer, allocatable, dimension(:), public :: nprkx
  integer, allocatable, dimension(:), public :: nhrkx
  integer, allocatable, dimension(:), public :: nlall
  integer, allocatable, dimension(:), public :: nlstst
  integer, allocatable, dimension(:), public :: nlf
  integer, allocatable, dimension(:), public :: nlmajmin
  integer, allocatable, dimension(:), public :: nldepd
  integer, allocatable, dimension(:), public :: nldepw
  integer, allocatable, dimension(:), public :: nldepx
  integer, allocatable, dimension(:,:), public :: njcoth
  integer, allocatable, dimension(:), public :: nmzjac
  integer, allocatable, dimension(:,:), public :: nzjac1
  integer, allocatable, dimension(:,:), public :: njcoss
  integer, allocatable, dimension(:), public :: nmsjac
  integer, allocatable, dimension(:,:), public :: nsjac1
  integer, allocatable, dimension(:), public :: nsspt
  integer, allocatable, dimension(:,:), public :: nspi
  integer, allocatable, dimension(:,:), public :: nsspi
  integer, allocatable, dimension(:), public :: nssi
  integer, allocatable, dimension(:), public :: nssrt
  integer, allocatable, dimension(:,:), public :: nssri
  integer, allocatable, dimension(:,:), public :: nssrx
  real, allocatable, dimension(:), public :: frpb
  real, allocatable, dimension(:), public :: frpt
  real, allocatable, dimension(:), public :: frpj
  real, allocatable, dimension(:), public :: frph
  real, allocatable, dimension(:), public :: frpx
  integer, allocatable, dimension(:,:), public :: nonzero_map_unordered
  integer, allocatable, dimension(:,:), public :: modified_map
  integer, allocatable, dimension(:,:), public :: nonzero_map
  integer, allocatable, dimension(:), public :: reorder
  integer, allocatable, dimension(:), public :: ilcf
  integer, allocatable, dimension(:), public :: ilss
  integer, allocatable, dimension(:), public :: ilct
  integer, allocatable, dimension(:), public :: ilftr
  integer, allocatable, dimension(:), public :: ilft
  integer, allocatable, dimension(:), public :: ilstmin
  logical, allocatable, dimension(:,:), public :: linfam
  logical, allocatable, dimension(:), public :: ldepd
  logical, allocatable, dimension(:), public :: ldepw
  character(len = 10), allocatable, public :: advt(:)
  character(len = 10), allocatable, public :: nadvt(:)
  character(len = 10), allocatable, public :: family(:)
  character(len = 10), allocatable, public :: speci(:)
  character(len = 10), allocatable, public :: spro2(:)
  character(len = 10), allocatable, public :: specf(:)
  character(len = 2), allocatable, public :: ctype(:)
  character(len = 10), allocatable, public :: spb(:, :)
  character(len = 10), allocatable, public :: spt(:, :)
  character(len = 10), allocatable, public :: spj(:, :)
  character(len = 10), allocatable, public :: sph(:, :)
  integer, public :: total
  integer, public :: jpfrpb
  integer, public :: jpfrpt
  integer, public :: jpfrpj
  integer, public :: jpfrph
  integer, public :: jpfrpx
  integer, save, public :: spfjsize_max
  integer, save, public :: maxterms
  integer, save, public :: maxfterms
  integer, allocatable, dimension(:), save, public :: nposterms
  integer, allocatable, dimension(:), save, public :: nnegterms
  integer, allocatable, dimension(:), save, public :: nfracterms
  integer, allocatable, dimension(:,:), save, public :: posterms
  integer, allocatable, dimension(:,:), save, public :: negterms
  integer, allocatable, dimension(:,:), save, public :: fracterms
  integer, allocatable, dimension(:), save, public :: base_tracer
  real, allocatable, dimension(:,:), save, public :: ffrac
  character(len = 2), parameter, public :: jpfm = 'FM'
  character(len = 2), parameter, public :: jpif = 'FT'
  character(len = 2), parameter, public :: jpsp = 'TR'
  character(len = 2), parameter, public :: jpna = 'SS'
  character(len = 2), parameter, public :: jpco = 'CT'
  character(len = 2), parameter, public :: jpcf = 'CF'
  character(len = 2), parameter, public :: jpoo = 'OO'
  logical, public :: o1d_in_ss
  logical, public :: o3p_in_ss
  logical, public :: n_in_ss
  logical, public :: h_in_ss
  real, public :: cdt
  real, public :: peps
  integer, public :: nitnr
  integer, public :: nitfg
  integer, public :: ntrf
  integer, public :: ntr3
  integer, public :: nnaf
  integer, public :: nro2
  integer, public :: nuni
  integer, public :: nsst
  integer, public :: ncsteps
  integer, save, public :: nit0 = 20
  integer, public :: nfphot
  integer, public :: jsubs
  integer, public :: method
  integer, save, public :: interval = imdi
  integer, public :: nnfrp
  integer, public :: nstst
  integer, public :: nf
  integer, public :: ndepd
  integer, public :: ndepw
  integer, public :: ntro3
  integer, public :: ntroh
  integer, public :: ntrho2
  integer, public :: ntrno
  integer, public :: nspo1d
  integer, public :: nspo3p
  integer, public :: nspo3
  integer, public :: nspoh
  integer, public :: nspho2
  integer, public :: nspno
  integer, public :: nspn
  integer, public :: nsph
  integer, public :: ih_o3
  integer, public :: ih_h2o2
  integer, public :: ih_so2
  integer, public :: ih_hno3
  integer, public :: ih_o3_const
  integer, public :: ihso3_h2o2
  integer, public :: ihso3_o3
  integer, public :: iso3_o3
  integer, public :: ih2so4_hv
  integer, public :: iso2_oh
  integer, public :: ih2o2_oh
  integer, public :: ihno3_oh
  integer, public :: in2o5_h
  integer, public :: iho2_h
  integer, public :: jsro2
  logical, public :: ljacx
  logical, public :: ltrig
  character(len = *), parameter, private :: ModuleName = 'ASAD_MOD'

  INTEGER, PARAMETER :: maxLineLen = 1024
  CHARACTER (LEN=maxLineLen)    :: umMessage

  public

  contains
  subroutine asad_mod_pre_setup_init()
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_MOD_PRE_SETUP_INIT'

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (.NOT.ALLOCATED(nodd)) then
      ALLOCATE(nodd(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nltrf)) then
      ALLOCATE(nltrf(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nltr3)) then
      ALLOCATE(nltr3(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nlnaro2)) then
      ALLOCATE(nlnaro2(1:jpro2))
    end if
    if (.NOT.ALLOCATED(nlfro2)) then
      ALLOCATE(nlfro2(1:jpro2))
    end if
    if (.NOT.ALLOCATED(advt)) then
      ALLOCATE(advt(1:jpctr))
    end if
    if (.NOT.ALLOCATED(nadvt)) then
      ALLOCATE(nadvt(1:jpspec - jpctr))
    end if
    if (.NOT.ALLOCATED(spro2)) then
      ALLOCATE(spro2(1:jpro2))
    end if
    if (.NOT.ALLOCATED(specf)) then
      ALLOCATE(specf(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(family)) then
      ALLOCATE(family(1:jpspec))
    end if
    if (.NOT.ALLOCATED(speci)) then
      ALLOCATE(speci(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ctype)) then
      ALLOCATE(ctype(1:jpspec))
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_mod_pre_setup_init
  subroutine asad_mod_init(n_points)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    logical, save :: firstcall = .true.
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_MOD_INIT'
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    if (firstcall) then
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Nullify_Stmt
      !  - Unsupported statement: Nullify_Stmt
      NULLIFY(prod)
      NULLIFY(slos)
      firstcall = .false.
    end if
    jpfrpb = (jpspb - 2) * jpbk
    jpfrpt = (jpspt - 2) * jptk
    jpfrpj = (jpspj - 2) * jppj
    jpfrph = (jpsph - 2) * jphk
    jpfrpx = jpfrpb + jpfrpt + jpfrpj + jpfrph
    if (.NOT.ALLOCATED(madvtr)) then
      ALLOCATE(madvtr(1:jpspec))
    end if
    if (.NOT.ALLOCATED(majors)) then
      ALLOCATE(majors(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(moffam)) then
      ALLOCATE(moffam(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ipa2)) then
      ALLOCATE(ipa2(1:jpctr))
    end if
    if (.NOT.ALLOCATED(nspi)) then
      ALLOCATE(nspi(1:jpnr,1:jpmsp))
    end if
    if (.NOT.ALLOCATED(nsspt)) then
      ALLOCATE(nsspt(1:jpss))
    end if
    if (.NOT.ALLOCATED(nsspi)) then
      ALLOCATE(nsspi(1:jpss,1:jpssr))
    end if
    if (.NOT.ALLOCATED(nssi)) then
      ALLOCATE(nssi(1:jpss))
    end if
    if (.NOT.ALLOCATED(nssrt)) then
      ALLOCATE(nssrt(1:jpss))
    end if
    if (.NOT.ALLOCATED(nssri)) then
      ALLOCATE(nssri(1:jpss,1:jpssr))
    end if
    if (.NOT.ALLOCATED(nssrx)) then
      ALLOCATE(nssrx(1:jpss,1:jpssr))
    end if
    if (.NOT.ALLOCATED(depvel)) then
      ALLOCATE(depvel(1:jddept,1:jddepc,1:jpdd))
    end if
    if (.NOT.ALLOCATED(k298)) then
      ALLOCATE(k298(1:jpdw))
    end if
    if (.NOT.ALLOCATED(dhr)) then
      ALLOCATE(dhr(1:jpdw))
    end if
    if (.NOT.ALLOCATED(kd298)) then
      ALLOCATE(kd298(1:jpdw,1:jpeq))
    end if
    if (.NOT.ALLOCATED(ddhr)) then
      ALLOCATE(ddhr(1:jpdw,1:jpeq))
    end if
    if (.NOT.ALLOCATED(ldepd)) then
      ALLOCATE(ldepd(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ldepw)) then
      ALLOCATE(ldepw(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nltrim)) then
      ALLOCATE(nltrim(0:jpcspf,1:3))
    end if
    if (.NOT.ALLOCATED(nlpdv)) then
      ALLOCATE(nlpdv(1:(jpspj - 2) * jppj,1:2))
    end if
    if (.NOT.ALLOCATED(ab)) then
      ALLOCATE(ab(1:jpbk,1:jpab))
    end if
    if (.NOT.ALLOCATED(at)) then
      ALLOCATE(at(1:jptk,1:jpat))
    end if
    if (.NOT.ALLOCATED(aj)) then
      ALLOCATE(aj(1:jppj,1:jpaj))
    end if
    if (.NOT.ALLOCATED(ah)) then
      ALLOCATE(ah(1:jphk,1:jpah))
    end if
    if (.NOT.ALLOCATED(spb)) then
      ALLOCATE(spb(1:jpbk + 1,1:jpspb))
    end if
    if (.NOT.ALLOCATED(spt)) then
      ALLOCATE(spt(1:jptk + 1,1:jpspt))
    end if
    if (.NOT.ALLOCATED(spj)) then
      ALLOCATE(spj(1:jppj + 1,1:jpspj))
    end if
    if (.NOT.ALLOCATED(sph)) then
      ALLOCATE(sph(1:jphk + 1,1:jpsph))
    end if
    if (.NOT.ALLOCATED(frpb)) then
      ALLOCATE(frpb(1:jpfrpb))
    end if
    if (.NOT.ALLOCATED(frpt)) then
      ALLOCATE(frpt(1:jpfrpt))
    end if
    if (.NOT.ALLOCATED(frpj)) then
      ALLOCATE(frpj(1:jpfrpj))
    end if
    if (.NOT.ALLOCATED(frph)) then
      ALLOCATE(frph(1:jpfrph))
    end if
    if (.NOT.ALLOCATED(frpx)) then
      ALLOCATE(frpx(1:jpfrpx))
    end if
    if (.NOT.ALLOCATED(ztabpd)) then
      ALLOCATE(ztabpd(1:jpfrpd,1:2))
    end if
    if (.NOT.ALLOCATED(nfrpx)) then
      ALLOCATE(nfrpx(1:jpnr))
    end if
    if (.NOT.ALLOCATED(ntabfp)) then
      ALLOCATE(ntabfp(1:jpfrpx,1:3))
    end if
    if (.NOT.ALLOCATED(ntabpd)) then
      ALLOCATE(ntabpd(1:jpfrpd,1:3))
    end if
    if (.NOT.ALLOCATED(npdfr)) then
      ALLOCATE(npdfr(1:jpnr,1:2))
    end if
    if (.NOT.ALLOCATED(ngrp)) then
      ALLOCATE(ngrp(1:2 * jpspec,1:3))
    end if
    if (.NOT.ALLOCATED(njcgrp)) then
      ALLOCATE(njcgrp(1:jpcspf,1:3))
    end if
    if (.NOT.ALLOCATED(nprdx3)) then
      ALLOCATE(nprdx3(1:3,1:jpnr / (3 * 3) + 3 * 3,1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(nprdx2)) then
      ALLOCATE(nprdx2(1:2,1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(nprdx1)) then
      ALLOCATE(nprdx1(1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(njacx3)) then
      ALLOCATE(njacx3(1:3,1:jpnr / (3 * 3) + 3 * 3,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(njacx2)) then
      ALLOCATE(njacx2(1:2,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(njacx1)) then
      ALLOCATE(njacx1(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nmpjac)) then
      ALLOCATE(nmpjac(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(npjac1)) then
      ALLOCATE(npjac1(1:jppjac,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nbrkx)) then
      ALLOCATE(nbrkx(1:jpbk + 1))
    end if
    if (.NOT.ALLOCATED(ntrkx)) then
      ALLOCATE(ntrkx(1:jptk + 1))
    end if
    if (.NOT.ALLOCATED(nprkx)) then
      ALLOCATE(nprkx(1:jppj + 1))
    end if
    if (.NOT.ALLOCATED(nhrkx)) then
      ALLOCATE(nhrkx(1:jphk + 1))
    end if
    if (.NOT.ALLOCATED(nlall)) then
      ALLOCATE(nlall(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nlstst)) then
      ALLOCATE(nlstst(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nlf)) then
      ALLOCATE(nlf(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nlmajmin)) then
      ALLOCATE(nlmajmin(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nldepd)) then
      ALLOCATE(nldepd(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nldepw)) then
      ALLOCATE(nldepw(1:jpspec))
    end if
    if (.NOT.ALLOCATED(nldepx)) then
      ALLOCATE(nldepx(1:jpspec + 6))
    end if
    if (.NOT.ALLOCATED(njcoth)) then
      ALLOCATE(njcoth(1:jpnr,1:jpmsp))
    end if
    if (.NOT.ALLOCATED(nmzjac)) then
      ALLOCATE(nmzjac(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nzjac1)) then
      ALLOCATE(nzjac1(1:jpnr,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(njcoss)) then
      ALLOCATE(njcoss(1:jpnr,1:jpmsp))
    end if
    if (.NOT.ALLOCATED(nmsjac)) then
      ALLOCATE(nmsjac(1:jpcspf))
    end if
    if (.NOT.ALLOCATED(nsjac1)) then
      ALLOCATE(nsjac1(1:jpnr,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(shno3)) then
      ALLOCATE(shno3(1:n_points))
    end if
    if (.NOT.ALLOCATED(sh2o)) then
      ALLOCATE(sh2o(1:n_points))
    end if
    if (.NOT.ALLOCATED(fpsc1)) then
      ALLOCATE(fpsc1(1:n_points))
    end if
    if (.NOT.ALLOCATED(fpsc2)) then
      ALLOCATE(fpsc2(1:n_points))
    end if
    if (.NOT.ALLOCATED(ilcf)) then
      ALLOCATE(ilcf(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ilss)) then
      ALLOCATE(ilss(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ilct)) then
      ALLOCATE(ilct(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ilftr)) then
      ALLOCATE(ilftr(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ilft)) then
      ALLOCATE(ilft(1:jpspec))
    end if
    if (.NOT.ALLOCATED(ilstmin)) then
      ALLOCATE(ilstmin(1:jpspec))
    end if
    !NO!$acc kernels
    method = ukca_config%ukca_int_method
    nitnr = 10
    nitfg = 10
    !NO!$acc loop gang vector independent collapse(2)
    do idx = LBOUND(njcoth, dim=2), UBOUND(njcoth, dim=2), 1
      do idx_1 = LBOUND(njcoth, dim=1), UBOUND(njcoth, dim=1), 1
        njcoth(idx_1,idx) = 0
      enddo
    enddo
    !NO!$acc end kernels
    if (method == int_method_nr) then
      if (.NOT.ALLOCATED(nonzero_map_unordered)) then
        ALLOCATE(nonzero_map_unordered(1:jpcspf,1:jpcspf))
      end if
      if (.NOT.ALLOCATED(modified_map)) then
        ALLOCATE(modified_map(1:jpcspf,1:jpcspf))
      end if
      if (.NOT.ALLOCATED(nonzero_map)) then
        ALLOCATE(nonzero_map(1:jpcspf,1:jpcspf))
      end if
      if (.NOT.ALLOCATED(reorder)) then
        ALLOCATE(reorder(1:jpcspf))
      end if
      if (ukca_config%l_ukca_cristrat) then
        spfjsize_max = 3000
        maxterms = 345
        maxfterms = 300
      else
        spfjsize_max = 1000
        maxterms = 160
        maxfterms = 100
      end if
      if (.NOT.ALLOCATED(nposterms)) then
        ALLOCATE(nposterms(1:spfjsize_max))
      end if
      if (.NOT.ALLOCATED(nnegterms)) then
        ALLOCATE(nnegterms(1:spfjsize_max))
      end if
      if (.NOT.ALLOCATED(nfracterms)) then
        ALLOCATE(nfracterms(1:spfjsize_max))
      end if
      if (.NOT.ALLOCATED(posterms)) then
        ALLOCATE(posterms(1:spfjsize_max,1:maxterms))
      end if
      if (.NOT.ALLOCATED(negterms)) then
        ALLOCATE(negterms(1:spfjsize_max,1:maxterms))
      end if
      if (.NOT.ALLOCATED(fracterms)) then
        ALLOCATE(fracterms(1:spfjsize_max,1:maxfterms))
      end if
      if (.NOT.ALLOCATED(base_tracer)) then
        ALLOCATE(base_tracer(1:spfjsize_max))
      end if
      if (.NOT.ALLOCATED(ffrac)) then
        ALLOCATE(ffrac(1:spfjsize_max,1:maxfterms))
      end if
    end if
    if (.NOT.ALLOCATED(pd)) then
      ALLOCATE(pd(1:n_points,1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(co3)) then
      ALLOCATE(co3(1:n_points))
    end if
    if (.NOT.ALLOCATED(deriv)) then
      ALLOCATE(deriv(1:n_points,1:4,1:4))
    end if
    if (.NOT.ALLOCATED(dpd)) then
      ALLOCATE(dpd(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(dpw)) then
      ALLOCATE(dpw(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ej)) then
      ALLOCATE(ej(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(emr)) then
      ALLOCATE(emr(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(f)) then
      ALLOCATE(f(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fdot)) then
      ALLOCATE(fdot(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fj)) then
      ALLOCATE(fj(1:n_points,1:jpcspf,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fpsc1)) then
      ALLOCATE(fpsc1(1:n_points))
    end if
    if (.NOT.ALLOCATED(fpsc2)) then
      ALLOCATE(fpsc2(1:n_points))
    end if
    if (.NOT.ALLOCATED(ftilde)) then
      ALLOCATE(ftilde(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(ipa)) then
      ALLOCATE(ipa(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(lati)) then
      ALLOCATE(lati(1:n_points))
    end if
    if (.NOT.ALLOCATED(linfam)) then
      ALLOCATE(linfam(1:n_points,0:jpcspf))
    end if
    if (.NOT.ALLOCATED(p)) then
      ALLOCATE(p(1:n_points))
    end if
    if (.NOT.ALLOCATED(pmintnd)) then
      ALLOCATE(pmintnd(1:n_points))
    end if
    if (.NOT.ALLOCATED(prk)) then
      ALLOCATE(prk(1:n_points,1:jpnr))
    end if
    if (.NOT.ALLOCATED(qa)) then
      ALLOCATE(qa(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ratio)) then
      ALLOCATE(ratio(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(rk)) then
      ALLOCATE(rk(1:n_points,1:jpnr))
    end if
    if (.NOT.ALLOCATED(sh2o)) then
      ALLOCATE(sh2o(1:n_points))
    end if
    if (.NOT.ALLOCATED(shno3)) then
      ALLOCATE(shno3(1:n_points))
    end if
    if (.NOT.ALLOCATED(sph2o)) then
      ALLOCATE(sph2o(1:n_points))
    end if
    if (.NOT.ALLOCATED(sphno3)) then
      ALLOCATE(sphno3(1:n_points))
    end if
    if (.NOT.ALLOCATED(t)) then
      ALLOCATE(t(1:n_points))
    end if
    if (.NOT.ALLOCATED(t300)) then
      ALLOCATE(t300(1:n_points))
    end if
    if (.NOT.ALLOCATED(tnd)) then
      ALLOCATE(tnd(1:n_points))
    end if
    if (.NOT.ALLOCATED(wp)) then
      ALLOCATE(wp(1:n_points))
    end if
    if (.NOT.ALLOCATED(co2)) then
      ALLOCATE(co2(1:n_points))
    end if
    if (.NOT.ALLOCATED(y)) then
      ALLOCATE(y(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ydot)) then
      ALLOCATE(ydot(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(za)) then
      ALLOCATE(za(1:n_points))
    end if
    if (method == int_method_nr) then
      if (.NOT.ALLOCATED(modified_map)) then
        ALLOCATE(modified_map(1:jpcspf,1:jpcspf))
      end if
      if (.NOT.ALLOCATED(spfj)) then
        ALLOCATE(spfj(1:n_points,1:spfjsize_max))
      end if
    end if

    prod => pd(:,1:jpspec)
    slos => pd(:,jpspec + 1:2 * jpspec)

    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do idx_2 = LBOUND(deriv, dim=3), UBOUND(deriv, dim=3), 1
      do idx_3 = LBOUND(deriv, dim=2), UBOUND(deriv, dim=2), 1
        do idx_4 = LBOUND(deriv, dim=1), UBOUND(deriv, dim=1), 1
          deriv(idx_4,idx_3,idx_2) = 1.0
        enddo
      enddo
    enddo
    !NO!$acc end kernels
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_mod_init
  subroutine asad_mod_init_spatial_vars(n_points)
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: n_points
    real(kind=jprb) :: zhook_handle
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_MOD_INIT_SPATIAL_VARS'
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

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    method = ukca_config%ukca_int_method
    if (ASSOCIATED(prod)) then
      print *, 'RW - nullifying prod'
      call umprint(ummessage, src='asad_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Nullify_Stmt
      NULLIFY(prod)
    end if
    if (ASSOCIATED(slos)) then
      print *, 'RW - nullifying slos'
      call umprint(ummessage, src='asad_mod')
      ! PSyclone CodeBlock (unsupported code) reason:
      !  - Unsupported statement: Nullify_Stmt
      NULLIFY(slos)
    end if
    if (.NOT.ALLOCATED(pd)) then
      print *, 'RW - allocating pd'
      call umprint(ummessage, src='asad_mod')
      ALLOCATE(pd(1:n_points,1:2 * jpspec))
    end if
    if (.NOT.ALLOCATED(co3)) then
      ALLOCATE(co3(1:n_points))
    end if
    if (.NOT.ALLOCATED(deriv)) then
      ALLOCATE(deriv(1:n_points,1:4,1:4))
    end if
    if (.NOT.ALLOCATED(dpd)) then
      ALLOCATE(dpd(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(dpw)) then
      ALLOCATE(dpw(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ej)) then
      ALLOCATE(ej(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(emr)) then
      ALLOCATE(emr(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(f)) then
      ALLOCATE(f(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fdot)) then
      ALLOCATE(fdot(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fj)) then
      ALLOCATE(fj(1:n_points,1:jpcspf,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(fpsc1)) then
      ALLOCATE(fpsc1(1:n_points))
    end if
    if (.NOT.ALLOCATED(fpsc2)) then
      ALLOCATE(fpsc2(1:n_points))
    end if
    if (.NOT.ALLOCATED(ftilde)) then
      ALLOCATE(ftilde(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(ipa)) then
      ALLOCATE(ipa(1:n_points,1:jpcspf))
    end if
    if (.NOT.ALLOCATED(lati)) then
      ALLOCATE(lati(1:n_points))
    end if
    if (.NOT.ALLOCATED(linfam)) then
      ALLOCATE(linfam(1:n_points,0:jpcspf))
    end if
    if (.NOT.ALLOCATED(p)) then
      ALLOCATE(p(1:n_points))
    end if
    if (.NOT.ALLOCATED(pmintnd)) then
      ALLOCATE(pmintnd(1:n_points))
    end if
    if (.NOT.ALLOCATED(prk)) then
      ALLOCATE(prk(1:n_points,1:jpnr))
    end if
    if (.NOT.ALLOCATED(qa)) then
      ALLOCATE(qa(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ratio)) then
      ALLOCATE(ratio(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(rk)) then
      ALLOCATE(rk(1:n_points,1:jpnr))
    end if
    if (.NOT.ALLOCATED(sh2o)) then
      ALLOCATE(sh2o(1:n_points))
    end if
    if (.NOT.ALLOCATED(shno3)) then
      ALLOCATE(shno3(1:n_points))
    end if
    if (.NOT.ALLOCATED(sph2o)) then
      ALLOCATE(sph2o(1:n_points))
    end if
    if (.NOT.ALLOCATED(sphno3)) then
      ALLOCATE(sphno3(1:n_points))
    end if
    if (.NOT.ALLOCATED(t)) then
      ALLOCATE(t(1:n_points))
    end if
    if (.NOT.ALLOCATED(t300)) then
      ALLOCATE(t300(1:n_points))
    end if
    if (.NOT.ALLOCATED(tnd)) then
      ALLOCATE(tnd(1:n_points))
    end if
    if (.NOT.ALLOCATED(wp)) then
      ALLOCATE(wp(1:n_points))
    end if
    if (.NOT.ALLOCATED(co2)) then
      ALLOCATE(co2(1:n_points))
    end if
    if (.NOT.ALLOCATED(y)) then
      ALLOCATE(y(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(ydot)) then
      ALLOCATE(ydot(1:n_points,1:jpspec))
    end if
    if (.NOT.ALLOCATED(za)) then
      ALLOCATE(za(1:n_points))
    end if

    print *, 'RW - LBOUND pd: ', LBOUND(pd, dim=1)
    call umprint(ummessage, src='asad_mod')
    print *, 'RW - UBOUND pd: ', ubound(pd, dim=1)
    call umprint(ummessage, src='asad_mod')

    if (method == int_method_nr) then
      if (.NOT.ALLOCATED(spfj)) then
        ALLOCATE(spfj(1:n_points,1:spfjsize_max))
      end if
      !NO!$acc kernels
      !NO!$acc loop gang vector independent collapse(2)
      do idx = LBOUND(spfj, dim=2), UBOUND(spfj, dim=2), 1
        do idx_1 = LBOUND(spfj, dim=1), UBOUND(spfj, dim=1), 1
          spfj(idx_1,idx) = 0.0
        enddo
      enddo
      !NO!$acc end kernels
    end if
    !NO!$acc kernels
    !NO!$acc loop gang vector independent collapse(3)
    do idx_2 = LBOUND(deriv, dim=3), UBOUND(deriv, dim=3), 1
      do idx_3 = LBOUND(deriv, dim=2), UBOUND(deriv, dim=2), 1
        do idx_4 = LBOUND(deriv, dim=1), UBOUND(deriv, dim=1), 1
          deriv(idx_4,idx_3,idx_2) = 1.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_5 = LBOUND(linfam, dim=2), UBOUND(linfam, dim=2), 1
      do idx_6 = LBOUND(linfam, dim=1), UBOUND(linfam, dim=1), 1
        linfam(idx_6,idx_5) = .false.
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_7 = LBOUND(co2, dim=1), UBOUND(co2, dim=1), 1
      co2(idx_7) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_8 = LBOUND(co3, dim=1), UBOUND(co3, dim=1), 1
      co3(idx_8) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_9 = LBOUND(dpd, dim=2), UBOUND(dpd, dim=2), 1
      do idx_10 = LBOUND(dpd, dim=1), UBOUND(dpd, dim=1), 1
        dpd(idx_10,idx_9) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_11 = LBOUND(dpw, dim=2), UBOUND(dpw, dim=2), 1
      do idx_12 = LBOUND(dpw, dim=1), UBOUND(dpw, dim=1), 1
        dpw(idx_12,idx_11) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_13 = LBOUND(ej, dim=2), UBOUND(ej, dim=2), 1
      do idx_14 = LBOUND(ej, dim=1), UBOUND(ej, dim=1), 1
        ej(idx_14,idx_13) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_15 = LBOUND(emr, dim=2), UBOUND(emr, dim=2), 1
      do idx_16 = LBOUND(emr, dim=1), UBOUND(emr, dim=1), 1
        emr(idx_16,idx_15) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_17 = LBOUND(f, dim=2), UBOUND(f, dim=2), 1
      do idx_18 = LBOUND(f, dim=1), UBOUND(f, dim=1), 1
        f(idx_18,idx_17) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_19 = LBOUND(fdot, dim=2), UBOUND(fdot, dim=2), 1
      do idx_20 = LBOUND(fdot, dim=1), UBOUND(fdot, dim=1), 1
        fdot(idx_20,idx_19) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(3)
    do idx_21 = LBOUND(fj, dim=3), UBOUND(fj, dim=3), 1
      do idx_22 = LBOUND(fj, dim=2), UBOUND(fj, dim=2), 1
        do idx_23 = LBOUND(fj, dim=1), UBOUND(fj, dim=1), 1
          fj(idx_23,idx_22,idx_21) = 0.0
        enddo
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_24 = LBOUND(fpsc1, dim=1), UBOUND(fpsc1, dim=1), 1
      fpsc1(idx_24) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_25 = LBOUND(fpsc2, dim=1), UBOUND(fpsc2, dim=1), 1
      fpsc2(idx_25) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_26 = LBOUND(ftilde, dim=2), UBOUND(ftilde, dim=2), 1
      do idx_27 = LBOUND(ftilde, dim=1), UBOUND(ftilde, dim=1), 1
        ftilde(idx_27,idx_26) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_28 = LBOUND(ipa, dim=2), UBOUND(ipa, dim=2), 1
      do idx_29 = LBOUND(ipa, dim=1), UBOUND(ipa, dim=1), 1
        ipa(idx_29,idx_28) = 0
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_30 = LBOUND(lati, dim=1), UBOUND(lati, dim=1), 1
      lati(idx_30) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_31 = LBOUND(p, dim=1), UBOUND(p, dim=1), 1
      p(idx_31) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_32 = LBOUND(pd, dim=2), UBOUND(pd, dim=2), 1
      do idx_33 = LBOUND(pd, dim=1), UBOUND(pd, dim=1), 1
        pd(idx_33,idx_32) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_34 = LBOUND(pmintnd, dim=1), UBOUND(pmintnd, dim=1), 1
      pmintnd(idx_34) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_37 = LBOUND(qa, dim=2), UBOUND(qa, dim=2), 1
      do idx_38 = LBOUND(qa, dim=1), UBOUND(qa, dim=1), 1
        qa(idx_38,idx_37) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_39 = LBOUND(ratio, dim=2), UBOUND(ratio, dim=2), 1
      do idx_40 = LBOUND(ratio, dim=1), UBOUND(ratio, dim=1), 1
        ratio(idx_40,idx_39) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_41 = LBOUND(sh2o, dim=1), UBOUND(sh2o, dim=1), 1
      sh2o(idx_41) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_42 = LBOUND(shno3, dim=1), UBOUND(shno3, dim=1), 1
      shno3(idx_42) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_45 = LBOUND(sph2o, dim=1), UBOUND(sph2o, dim=1), 1
      sph2o(idx_45) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_46 = LBOUND(sphno3, dim=1), UBOUND(sphno3, dim=1), 1
      sphno3(idx_46) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_47 = LBOUND(t, dim=1), UBOUND(t, dim=1), 1
      t(idx_47) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_48 = LBOUND(t300, dim=1), UBOUND(t300, dim=1), 1
      t300(idx_48) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_49 = LBOUND(tnd, dim=1), UBOUND(tnd, dim=1), 1
      tnd(idx_49) = 0.0
    enddo
    !NO!$acc loop gang vector independent
    do idx_50 = LBOUND(wp, dim=1), UBOUND(wp, dim=1), 1
      wp(idx_50) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_51 = LBOUND(y, dim=2), UBOUND(y, dim=2), 1
      do idx_52 = LBOUND(y, dim=1), UBOUND(y, dim=1), 1
        y(idx_52,idx_51) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_53 = LBOUND(ydot, dim=2), UBOUND(ydot, dim=2), 1
      do idx_54 = LBOUND(ydot, dim=1), UBOUND(ydot, dim=1), 1
        ydot(idx_54,idx_53) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent
    do idx_55 = LBOUND(za, dim=1), UBOUND(za, dim=1), 1
      za(idx_55) = 0.0
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_56 = LBOUND(rk, dim=2), UBOUND(rk, dim=2), 1
      do idx_57 = LBOUND(rk, dim=1), UBOUND(rk, dim=1), 1
        rk(idx_57,idx_56) = 0.0
      enddo
    enddo
    !NO!$acc loop gang vector independent collapse(2)
    do idx_58 = LBOUND(prk, dim=2), UBOUND(prk, dim=2), 1
      do idx_59 = LBOUND(prk, dim=1), UBOUND(prk, dim=1), 1
        prk(idx_59,idx_58) = 0.0
      enddo
    enddo
    !NO!$acc end kernels

    prod => pd(:,1:jpspec)
    slos => pd(:,jpspec + 1:2 * jpspec)

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()
    return

  end subroutine asad_mod_init_spatial_vars
  subroutine asad_mod_dealloc_spatial_vars()
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'ASAD_MOD_DEALLOC_SPATIAL_VARS'
    real(kind=jprb) :: zhook_handle

    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_in, zhook_handle)
    end if
    call nvtxStartRange(ModuleName // ':' // RoutineName)
    method = ukca_config%ukca_int_method
    if (method == int_method_nr) then
      if (ALLOCATED(spfj)) then
        DEALLOCATE(spfj)
      end if
    end if
    if (ALLOCATED(za)) then
      DEALLOCATE(za)
    end if
    if (ALLOCATED(ydot)) then
      DEALLOCATE(ydot)
    end if
    if (ALLOCATED(y)) then
      DEALLOCATE(y)
    end if
    if (ALLOCATED(co2)) then
      DEALLOCATE(co2)
    end if
    if (ALLOCATED(wp)) then
      DEALLOCATE(wp)
    end if
    if (ALLOCATED(tnd)) then
      DEALLOCATE(tnd)
    end if
    if (ALLOCATED(t300)) then
      DEALLOCATE(t300)
    end if
    if (ALLOCATED(t)) then
      DEALLOCATE(t)
    end if
    if (ALLOCATED(sphno3)) then
      DEALLOCATE(sphno3)
    end if
    if (ALLOCATED(sph2o)) then
      DEALLOCATE(sph2o)
    end if
    if (ALLOCATED(shno3)) then
      DEALLOCATE(shno3)
    end if
    if (ALLOCATED(sh2o)) then
      DEALLOCATE(sh2o)
    end if
    if (ALLOCATED(rk)) then
      DEALLOCATE(rk)
    end if
    if (ALLOCATED(ratio)) then
      DEALLOCATE(ratio)
    end if
    if (ALLOCATED(qa)) then
      DEALLOCATE(qa)
    end if
    if (ALLOCATED(prk)) then
      DEALLOCATE(prk)
    end if
    if (ALLOCATED(pmintnd)) then
      DEALLOCATE(pmintnd)
    end if
    if (ALLOCATED(p)) then
      DEALLOCATE(p)
    end if
    if (ALLOCATED(linfam)) then
      DEALLOCATE(linfam)
    end if
    if (ALLOCATED(lati)) then
      DEALLOCATE(lati)
    end if
    if (ALLOCATED(ipa)) then
      DEALLOCATE(ipa)
    end if
    if (ALLOCATED(ftilde)) then
      DEALLOCATE(ftilde)
    end if
    if (ALLOCATED(fpsc2)) then
      DEALLOCATE(fpsc2)
    end if
    if (ALLOCATED(fpsc1)) then
      DEALLOCATE(fpsc1)
    end if
    if (ALLOCATED(fj)) then
      DEALLOCATE(fj)
    end if
    if (ALLOCATED(fdot)) then
      DEALLOCATE(fdot)
    end if
    if (ALLOCATED(f)) then
      DEALLOCATE(f)
    end if
    if (ALLOCATED(emr)) then
      DEALLOCATE(emr)
    end if
    if (ALLOCATED(ej)) then
      DEALLOCATE(ej)
    end if
    if (ALLOCATED(dpw)) then
      DEALLOCATE(dpw)
    end if
    if (ALLOCATED(dpd)) then
      DEALLOCATE(dpd)
    end if
    if (ALLOCATED(deriv)) then
      DEALLOCATE(deriv)
    end if
    if (ALLOCATED(co3)) then
      DEALLOCATE(co3)
    end if
    if (ALLOCATED(pd)) then
      print *, 'RW - deallocating pd'
      call umprint(ummessage, src='asad_mod')
      DEALLOCATE(pd)
    end if
    if (lhook) then
      call dr_hook(ModuleName // ':' // RoutineName, zhook_out, zhook_handle)
    end if
    call nvtxEndRange()

  end subroutine asad_mod_dealloc_spatial_vars

end module asad_mod
