!    -------------------------------------------------------------------
!    Subroutine: READ_ENVIRONMENT_PARS_MOD                --------------
!
!    Purpose: Programme to read environmental parameters namelist
!             for UKCA box model.
!
!    Programming standard: @£($@*£$)
!
!    External documentation: @£*£@£%
!
!    -------------------------------------------------------------------
!
! Module for UKCA Box model namelist : ENVIRONMENT_PARS

MODULE read_environment_pars_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_ENVIRONMENT_PARS'

! Match to equivalent variables in code without "_1"
! Note names of variables must be different for some compilers to work
LOGICAL :: land_sea_mask_1
REAL    :: latitude_1
REAL    :: longitude_1
REAL    :: conv_cloud_lwp_1
REAL    :: tstar_1
REAL    :: zbl_1
REAL    :: rough_length_1
REAL    :: seaice_frac_1
REAL    :: pstar_1
REAL    :: zhsc_1
REAL    :: u_scalar_10m_1
REAL    :: u_s_1
REAL    :: surf_albedo_1
REAL    :: dms_sea_conc_1
INTEGER :: kent_1
INTEGER :: kent_dsc_1
INTEGER :: conv_cloud_base_1
INTEGER :: conv_cloud_top_1
REAL    :: q_1
REAL    :: qcf_1
REAL    :: qcl_1
REAL    :: conv_cloud_amount_1
REAL    :: rho_r2_1
REAL    :: area_cloud_fraction_1
REAL    :: cloud_frac_1
REAL    :: cloud_liq_frac_1
REAL    :: p_theta_levels_1
REAL    :: t_theta_levels_1
REAL    :: rhokh_mix_1
REAL    :: dtrdz_charney_grid_1
REAL    :: rhokh_rdz_1
REAL    :: dtrdz_1
REAL    :: we_lim_1
REAL    :: t_frac_1
REAL    :: zrzi_1
REAL    :: we_lim_dsc_1
REAL    :: t_frac_dsc_1
REAL    :: zrzi_dsc_1
REAL    :: ls_rain3d_1
REAL    :: ls_snow3d_1
REAL    :: autoconv_1
REAL    :: accretion_1
REAL    :: pv_on_theta_mlevs_1
REAL    :: conv_rain3d_1
REAL    :: conv_snow3d_1
REAL    :: so4_sa_clim_1

! Public procedures
PUBLIC :: read_environment_pars

!---------------------------------------------------------------------------
! Define namelist
!---------------------------------------------------------------------------

NAMELIST/ENVIRONMENT_PARS/ latitude_1, longitude_1, land_sea_mask_1,           &
  conv_cloud_lwp_1, tstar_1, zbl_1, rough_length_1, seaice_frac_1,             &
  pstar_1, zhsc_1, u_scalar_10m_1, u_s_1, surf_albedo_1, dms_sea_conc_1,       &
  kent_1, kent_dsc_1, conv_cloud_base_1, conv_cloud_top_1, q_1, qcf_1, qcl_1,  &
  conv_cloud_amount_1, rho_r2_1, area_cloud_fraction_1, cloud_frac_1,          &
  cloud_liq_frac_1, p_theta_levels_1, t_theta_levels_1, rhokh_mix_1,           &
  dtrdz_charney_grid_1, rhokh_rdz_1, dtrdz_1, we_lim_1, t_frac_1, zrzi_1,      &
  we_lim_dsc_1, t_frac_dsc_1, zrzi_dsc_1, ls_rain3d_1, ls_snow3d_1, autoconv_1,&
  accretion_1, pv_on_theta_mlevs_1, conv_rain3d_1, conv_snow3d_1,              &
  so4_sa_clim_1

PRIVATE :: ENVIRONMENT_PARS


! ----------------------------------------------------------------------

CONTAINS

SUBROUTINE read_environment_pars(icode, iomessage, ukca_box_nml,              &
  latitude_in, longitude_in, land_sea_mask_in, conv_cloud_lwp_in, tstar_in,   &
  zbl_in, rough_length_in, seaice_frac_in, pstar_in, zhsc_in,                 &
  u_scalar_10m_in, u_s_in, surf_albedo_in, dms_sea_conc_in, kent_in,          &
  kent_dsc_in, conv_cloud_base_in, conv_cloud_top_in, q_in, qcf_in, qcl_in,   &
  conv_cloud_amount_in, rho_r2_in, area_cloud_fraction_in, cloud_frac_in,     &
  cloud_liq_frac_in, p_theta_levels_in, t_theta_levels_in, rhokh_mix_in,      &
  dtrdz_charney_grid_in, rhokh_rdz_in, dtrdz_in, we_lim_in, t_frac_in,        &
  zrzi_in, we_lim_dsc_in, t_frac_dsc_in, zrzi_dsc_in, ls_rain3d_in,           &
  ls_snow3d_in, autoconv_in, accretion_in, pv_on_theta_mlevs_in,              &
  conv_rain3d_in, conv_snow3d_in, so4_sa_clim_in)

USE filenamelength_mod, ONLY: filenamelength

USE missing_data_mod,       ONLY: rmdi, imdi

USE umPrintMgr,             ONLY: umMessage, umPrint, PrintStatus, PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

!=============================================================================

IMPLICIT NONE

INTEGER, INTENT(OUT)                           :: icode
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: iomessage
CHARACTER(LEN=filenamelength), INTENT(IN)      :: ukca_box_nml

!++SAN scalar versions of variables for reading in namelist input
REAL, INTENT(OUT)    :: latitude_in
REAL, INTENT(OUT)    :: longitude_in
LOGICAL, INTENT(OUT) :: land_sea_mask_in
REAL, INTENT(OUT)    :: conv_cloud_lwp_in
REAL, INTENT(OUT)    :: tstar_in
REAL, INTENT(OUT)    :: zbl_in
REAL, INTENT(OUT)    :: rough_length_in
REAL, INTENT(OUT)    :: seaice_frac_in
REAL, INTENT(OUT)    :: pstar_in
REAL, INTENT(OUT)    :: zhsc_in
REAL, INTENT(OUT)    :: u_scalar_10m_in
REAL, INTENT(OUT)    :: u_s_in
REAL, INTENT(OUT)    :: surf_albedo_in
REAL, INTENT(OUT)    :: dms_sea_conc_in
INTEGER, INTENT(OUT) :: kent_in
INTEGER, INTENT(OUT) :: kent_dsc_in
INTEGER, INTENT(OUT) :: conv_cloud_base_in
INTEGER, INTENT(OUT) :: conv_cloud_top_in
REAL, INTENT(OUT)    :: q_in
REAL, INTENT(OUT)    :: qcf_in
REAL, INTENT(OUT)    :: qcl_in
REAL, INTENT(OUT)    :: conv_cloud_amount_in
REAL, INTENT(OUT)    :: rho_r2_in
REAL, INTENT(OUT)    :: area_cloud_fraction_in
REAL, INTENT(OUT)    :: cloud_frac_in
REAL, INTENT(OUT)    :: cloud_liq_frac_in
REAL, INTENT(OUT)    :: p_theta_levels_in
REAL, INTENT(OUT)    :: t_theta_levels_in
REAL, INTENT(OUT)    :: rhokh_mix_in
REAL, INTENT(OUT)    :: dtrdz_charney_grid_in
REAL, INTENT(OUT)    :: rhokh_rdz_in
REAL, INTENT(OUT)    :: dtrdz_in
REAL, INTENT(OUT)    :: we_lim_in
REAL, INTENT(OUT)    :: t_frac_in
REAL, INTENT(OUT)    :: zrzi_in
REAL, INTENT(OUT)    :: we_lim_dsc_in
REAL, INTENT(OUT)    :: t_frac_dsc_in
REAL, INTENT(OUT)    :: zrzi_dsc_in
REAL, INTENT(OUT)    :: ls_rain3d_in
REAL, INTENT(OUT)    :: ls_snow3d_in
REAL, INTENT(OUT)    :: autoconv_in
REAL, INTENT(OUT)    :: accretion_in
REAL, INTENT(OUT)    :: pv_on_theta_mlevs_in
REAL, INTENT(OUT)    :: conv_rain3d_in
REAL, INTENT(OUT)    :: conv_snow3d_in
REAL, INTENT(OUT)    :: so4_sa_clim_in

! Local variables
INTEGER :: istatus
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_ENVIROMENT_PARS'

!=============================================================================

! Initialise values
latitude_1              = rmdi
longitude_1             = rmdi
land_sea_mask_1         = .FALSE.
conv_cloud_lwp_1        = rmdi
tstar_1                 = rmdi
zbl_1                   = rmdi
rough_length_1          = rmdi
seaice_frac_1           = rmdi
pstar_1                 = rmdi
zhsc_1                  = rmdi
u_scalar_10m_1          = rmdi
u_s_1                   = rmdi
surf_albedo_1           = rmdi
dms_sea_conc_1          = rmdi
kent_1                  = imdi
kent_dsc_1              = imdi
conv_cloud_base_1       = imdi
conv_cloud_top_1        = imdi
q_1                     = rmdi
qcf_1                   = rmdi
qcl_1                   = rmdi
conv_cloud_amount_1     = rmdi
rho_r2_1                = rmdi
area_cloud_fraction_1   = rmdi
cloud_frac_1            = rmdi
cloud_liq_frac_1        = rmdi
p_theta_levels_1        = rmdi
t_theta_levels_1        = rmdi
rhokh_mix_1             = rmdi
dtrdz_charney_grid_1    = rmdi
rhokh_rdz_1             = rmdi
dtrdz_1                 = rmdi
we_lim_1                = rmdi
t_frac_1                = rmdi
zrzi_1                  = rmdi
we_lim_dsc_1            = rmdi
t_frac_dsc_1            = rmdi
zrzi_dsc_1              = rmdi
ls_rain3d_1             = rmdi
ls_snow3d_1             = rmdi
autoconv_1              = rmdi
accretion_1             = rmdi
pv_on_theta_mlevs_1     = rmdi
conv_rain3d_1           = rmdi
conv_snow3d_1           = rmdi
so4_sa_clim_1           = rmdi

!++SAN Setting up to read from a namelist
IF (PrintStatus >= PrStatus_Oper) THEN
  WRITE(umMessage,*) 'SAN; READ_ENVIRONMENT_PARS: Opening box model namelist'
  CALL umPrint(umMessage, src='read_environment_pars')
  WRITE(umMessage,*) 'Namelist file called ', ukca_box_nml
  CALL umPrint(umMessage, src='read_environment_pars')
END IF

OPEN(67,FILE=TRIM(ADJUSTL(ukca_box_nml)), ACTION='READ', IOSTAT=istatus)

IF (istatus /= 0) THEN
  icode = 500
  iomessage = ' Error opening file on unit 67 from '//routinename
  WRITE(umMessage,*) iomessage
  CALL umPrint(umMessage,src='read_environment_pars')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(ukca_box_nml))
  CALL umPrint(umMessage,src='read_environment_pars')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='read_environment_pars')
ELSE
  !++PMJ check initial values of variable
  WRITE(umMessage,*) 'PMJ; READING NAMELIST '
  CALL umPrint(umMessage, src='read_environment_pars')
  !++PMJ check  initial values of variable 
  READ (UNIT=67, NML=ENVIRONMENT_PARS)
END IF

CLOSE (67)

! Map variables
latitude_in            = latitude_1
longitude_in           = longitude_1
land_sea_mask_in       = land_sea_mask_1
conv_cloud_lwp_in      = conv_cloud_lwp_1
tstar_in               = tstar_1
zbl_in                 = zbl_1
rough_length_in        = rough_length_1
seaice_frac_in         = seaice_frac_1
pstar_in               = pstar_1
zhsc_in                = zhsc_1
u_scalar_10m_in        = u_scalar_10m_1
u_s_in                 = u_s_1
surf_albedo_in         = surf_albedo_1
kent_in                = kent_1
kent_dsc_in            = kent_dsc_1
q_in                   = q_1
qcf_in                 = qcf_1
qcl_in                 = qcl_1
conv_cloud_amount_in   = conv_cloud_amount_1
rho_r2_in              = rho_r2_1
area_cloud_fraction_in = area_cloud_fraction_1
cloud_frac_in          = cloud_frac_1
cloud_liq_frac_in      = cloud_liq_frac_1
p_theta_levels_in      = p_theta_levels_1
t_theta_levels_in      = t_theta_levels_1
rhokh_mix_in           = rhokh_mix_1
dtrdz_charney_grid_in  = dtrdz_charney_grid_1
rhokh_rdz_in           = rhokh_rdz_1
dtrdz_in               = dtrdz_1
we_lim_in              = we_lim_1
t_frac_in              = t_frac_1
zrzi_in                = zrzi_1
we_lim_dsc_in          = we_lim_dsc_1
t_frac_dsc_in          = t_frac_dsc_1
zrzi_dsc_in            = zrzi_dsc_1
ls_rain3d_in           = ls_rain3d_1
ls_snow3d_in           = ls_snow3d_1
autoconv_in            = autoconv_1
accretion_in           = accretion_1
pv_on_theta_mlevs_in   = pv_on_theta_mlevs_1
conv_cloud_base_in     = conv_cloud_base_1
conv_cloud_top_in      = conv_cloud_top_1
conv_rain3d_in         = conv_rain3d_1
conv_snow3d_in         = conv_snow3d_1
so4_sa_clim_in         = so4_sa_clim_1
dms_sea_conc_in        = dms_sea_conc_1

!++PMJ check initial values of variable
WRITE(umMessage,*) 'PMJ; READ_ENVIRONMENT_PARS: land_sea_mask = ', land_sea_mask_in
CALL umPrint(umMessage, src='read_environment_pars')
WRITE(umMessage,*) 'PMJ; READ_ENVIRONMENT_PARS: t_theta_levels = ', t_theta_levels_in
CALL umPrint(umMessage, src='read_environment_pars')
WRITE(umMessage,*) 'PMJ; READ_ENVIRONMENT_PARS: longitude and latitude= ', longitude_in, latitude_in
CALL umPrint(umMessage, src='read_environment_pars')
!++PMJ check  initial values of variable


RETURN
END SUBROUTINE read_environment_pars

!=============================================================================
END MODULE read_environment_pars_mod

