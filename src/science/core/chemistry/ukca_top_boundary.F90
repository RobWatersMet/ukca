! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Routine to impose a top boundary condition for species that have
!  a thermospheric source or sink, i.e. NO, CO, H2O.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_topboundary_mod

   IMPLICIT NONE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_TOPBOUNDARY_MOD'

CONTAINS

   SUBROUTINE ukca_top_boundary(row_length, rows, model_levels, ntracers, &
                                latitude, tracer)

      USE ukca_config_specification_mod, ONLY: ukca_config, i_top_BC_H2O
      USE asad_mod, ONLY: advt, jpctr, peps
      USE ukca_constants, ONLY: pi, c_no, c_co, c_h2o, c_h2, c_h, c_oh, c_o3
      USE ukca_time_mod, ONLY: i_day_number, i_year, days_in_year

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: ntracers
      REAL, INTENT(IN) :: latitude(:)  ! Latitude of each row (degrees)
      REAL, INTENT(IN OUT) :: tracer(row_length, rows, model_levels, ntracers)

! local variables
! Positions of NO and CO in UKCA tracers array
      INTEGER, SAVE :: i_no = -1
      INTEGER, SAVE :: i_co = -1
      INTEGER, SAVE :: i_h2o = -1
      INTEGER, SAVE :: i_h2 = -1
      INTEGER, SAVE :: i_h = -1
      INTEGER, SAVE :: i_oh = -1
      INTEGER, SAVE :: i_o3 = -1

      INTEGER, PARAMETER :: n_top = 4
      INTEGER, PARAMETER :: j_no = 1
      INTEGER, PARAMETER :: j_co = 2
      INTEGER, PARAMETER :: j_h2o = 3
      INTEGER, PARAMETER :: j_o3 = 4
! Number of latitudes in ACE FTS dataset
      INTEGER, PARAMETER :: n_ace = 36
! Latitude spacing in ACE-FTS dataset
      REAL, PARAMETER :: dellat_ace = 180.0/n_ace

      LOGICAL, SAVE :: firstcall = .TRUE.
      INTEGER, SAVE :: year_filled = 0

      REAL, ALLOCATABLE :: clim(:, :, :)
      REAL, ALLOCATABLE, SAVE :: clim_interp(:, :, :)
      REAL :: acelat(n_ace)
      REAL :: frac
      INTEGER :: i
      INTEGER :: j
      INTEGER :: idx
      INTEGER :: yearlen
      REAL, ALLOCATABLE :: change_h2o(:)
      REAL, ALLOCATABLE :: toth(:)

!ACE-FTS 3-monthly CO dataset given as constant and 1st and 2nd.
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Coeffs 1 2 3 4 5
      REAL, PARAMETER :: cocoeffs(5, n_ace) = RESHAPE([ &
                                                      1.0343E-05, -5.7139E-06, -1.0139E-06, -3.0417E-06, -1.5135E-06, &
                                                      1.0351E-05, -5.7028E-06, -1.0336E-06, -3.0498E-06, -1.5180E-06, &
                                                      9.7896E-06, -5.6210E-06, -1.3558E-06, -2.4121E-06, -1.1959E-06, &
                                                      9.3485E-06, -5.5150E-06, -2.1491E-06, -1.7859E-06, -8.8550E-07, &
                                                      6.4714E-06, -5.7903E-06, 1.6208E-06, 6.8119E-07, 4.2493E-07, &
                                                      6.5282E-06, -5.2432E-06, 9.1142E-07, 5.6361E-07, 3.5104E-07, &
                                                      6.7975E-06, -4.6991E-06, 5.3621E-07, 1.8180E-07, 1.4572E-07, &
                                                      7.3118E-06, -4.4828E-06, -2.9673E-08, 2.4394E-07, 1.6882E-07, &
                                                      7.6905E-06, -4.0141E-06, -4.2918E-07, 1.7204E-07, 1.2271E-07, &
                                                      8.3275E-06, -3.7961E-06, -7.4087E-07, -1.6181E-07, -5.4085E-08, &
                                                      8.7063E-06, -2.9312E-06, -1.3273E-06, -5.7266E-07, -2.7991E-07, &
                                                      9.6419E-06, -2.2294E-06, -9.7496E-07, -1.6869E-07, -7.5792E-08, &
                                                      9.8245E-06, -1.7929E-06, -8.1412E-07, -8.7181E-07, -4.3839E-07, &
                                                      1.0560E-05, -1.1445E-06, 1.0073E-07, -2.1751E-06, -1.1017E-06, &
                                                      1.1308E-05, -8.0858E-07, 4.9043E-07, -2.2772E-06, -1.1528E-06, &
                                                      1.1124E-05, -3.1799E-07, 1.4939E-07, -2.3904E-06, -1.2196E-06, &
                                                      1.1294E-05, 7.5480E-08, 2.2989E-08, -2.1227E-06, -1.0878E-06, &
                                                      1.1266E-05, 5.2460E-09, 5.1398E-07, -2.2039E-06, -1.1230E-06, &
                                                      1.1335E-05, -2.2720E-07, -6.5253E-08, -2.6727E-06, -1.3676E-06, &
                                                      1.0933E-05, -4.5845E-08, -1.8094E-07, -2.6850E-06, -1.3770E-06, &
                                                      1.0745E-05, 1.0299E-07, 2.8773E-08, -2.7200E-06, -1.3940E-06, &
                                                      1.0618E-05, 6.9697E-07, -4.1458E-07, -2.3916E-06, -1.2368E-06, &
                                                      9.9231E-06, 1.2579E-06, -1.9564E-08, -2.1797E-06, -1.1292E-06, &
                                                      9.1729E-06, 2.0863E-06, 9.3015E-07, -9.2643E-07, -4.8428E-07, &
                                                      8.3174E-06, 1.9475E-06, 1.5553E-06, -3.5080E-07, -1.8076E-07, &
                                                      7.9472E-06, 2.9051E-06, 1.6206E-06, 4.3667E-07, 2.1392E-07, &
                                                      8.0176E-06, 3.8365E-06, 1.9448E-06, 4.5154E-07, 2.1614E-07, &
                                                      7.4619E-06, 4.5268E-06, 1.2280E-06, 7.5508E-07, 3.5646E-07, &
                                                      6.9411E-06, 4.8276E-06, 7.1664E-07, 9.7551E-07, 4.6044E-07, &
                                                      6.6307E-06, 5.5499E-06, -3.9214E-08, 1.3711E-06, 6.4714E-07, &
                                                      6.1959E-06, 5.5581E-06, -3.2250E-07, 1.4772E-06, 6.9812E-07, &
                                                      6.1178E-06, 6.0612E-06, -1.0747E-06, 1.4358E-06, 6.6319E-07, &
                                                      9.2436E-06, 5.1172E-06, 3.9615E-06, -1.8447E-06, -9.4921E-07, &
                                                      9.4305E-06, 5.1483E-06, 3.3518E-06, -2.2105E-06, -1.1440E-06, &
                                                      1.0374E-05, 5.4509E-06, 2.1734E-06, -3.2752E-06, -1.7061E-06, &
                                                      1.0585E-05, 5.5457E-06, 1.8000E-06, -3.5124E-06, -1.8329E-06], &
                                                      [5, n_ace])

!ACE-FTS 3-monthly NO dataset decomposed into constant, 1st and 2nd
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Latitude Coeffs 1 2 3 4
      REAL, PARAMETER :: nocoeffs(5, n_ace) = RESHAPE([ &
                                                      3.7216E-07, -4.3727E-07, -5.5703E-08, 1.6251E-07, 8.6909E-08, &
                                                      3.5485E-07, -4.2486E-07, -7.9603E-08, 1.7401E-07, 9.2396E-08, &
                                                      3.3461E-07, -4.0067E-07, -1.2452E-07, 1.8099E-07, 9.5208E-08, &
                                                      2.8363E-07, -3.4596E-07, -8.8981E-08, 1.6209E-07, 8.5403E-08, &
                                                      4.7979E-07, -5.4090E-07, 1.7642E-07, 1.0260E-07, 5.9949E-08, &
                                                      4.4516E-07, -5.0789E-07, 1.1355E-07, 1.2274E-07, 6.9204E-08, &
                                                      2.5505E-07, -2.9443E-07, 3.2521E-08, 9.6749E-08, 5.2841E-08, &
                                                      1.2788E-07, -1.5069E-07, -2.8978E-08, 8.1171E-08, 4.2727E-08, &
                                                      9.2971E-08, -9.9975E-08, -1.7173E-08, 5.5500E-08, 2.9215E-08, &
                                                      9.8711E-08, -1.0974E-07, -2.0047E-08, 5.9627E-08, 3.1392E-08, &
                                                      5.7713E-08, -3.1688E-08, 1.2337E-09, 1.6967E-08, 9.0179E-09, &
                                                      4.7535E-08, -1.4523E-08, 4.5995E-10, 8.7556E-09, 4.6335E-09, &
                                                      4.6706E-08, -1.3386E-08, 7.6621E-09, 3.5865E-09, 2.0585E-09, &
                                                      4.1139E-08, -4.4998E-09, -2.0247E-09, -2.0725E-09, -1.0409E-09, &
                                                      4.1288E-08, -1.5117E-08, 3.0356E-09, -4.9061E-09, -2.3286E-09, &
                                                      4.4847E-08, 2.8043E-09, 1.0589E-08, -2.8860E-09, -1.3825E-09, &
                                                      4.3541E-08, -8.4126E-09, 4.8296E-10, 3.8074E-09, 2.0389E-09, &
                                                      4.7605E-08, -5.2713E-09, -1.3116E-09, 3.4990E-09, 1.8290E-09, &
                                                      5.1771E-08, 7.5440E-09, -7.9053E-09, -3.5942E-09, -2.0077E-09, &
                                                      4.1352E-08, -6.0758E-10, -7.5309E-09, -7.3710E-09, -3.8575E-09, &
                                                      4.3036E-08, 4.1618E-09, -1.1526E-08, -5.0942E-09, -2.7849E-09, &
                                                      4.1913E-08, 9.4329E-09, -1.9973E-08, -1.0213E-08, -5.5573E-09, &
                                                      4.0121E-08, 9.3359E-09, -8.6942E-09, -1.9265E-09, -1.1803E-09, &
                                                      4.5106E-08, 1.5444E-08, -4.6932E-09, 3.6038E-09, 1.6390E-09, &
                                                      4.6951E-08, 2.7356E-08, 5.5145E-10, 1.2740E-08, 6.2626E-09, &
                                                      4.7500E-08, 2.9097E-08, 5.8535E-09, 1.8281E-08, 9.1457E-09, &
                                                      6.2896E-08, 5.0146E-08, 8.3025E-09, 2.9605E-08, 1.4767E-08, &
                                                      8.2373E-08, 7.1688E-08, 2.9553E-09, 3.7750E-08, 1.8665E-08, &
                                                      1.0801E-07, 1.1840E-07, 1.2458E-08, 6.7333E-08, 3.3468E-08, &
                                                      1.8258E-07, 2.0183E-07, -2.3098E-08, 7.8927E-08, 3.8169E-08, &
                                                      2.5703E-07, 2.9406E-07, -5.0007E-08, 9.5218E-08, 4.5291E-08, &
                                                      2.2001E-07, 2.2272E-07, -5.9467E-08, 5.5834E-08, 2.5711E-08, &
                                                      2.1587E-07, 2.1169E-07, -9.0752E-09, 6.1104E-08, 2.9106E-08, &
                                                      2.3610E-07, 2.2014E-07, 5.0957E-08, 6.5617E-08, 3.2035E-08, &
                                                      1.7327E-07, 1.1810E-07, 5.4228E-08, 1.1941E-08, 5.5842E-09, &
                                                      1.9367E-07, 1.2734E-07, 1.8236E-08, -1.0867E-08, -6.6098E-09], &
                                                      [5, n_ace])

!ACE-FTS 3-monthly H2O dataset decomposed into constant, 1st and 2nd
!Fourier transform coefficients.
!Olaf Morgenstern, Sept 2018
!Latitude Coeffs 1 2 3 4
      REAL, PARAMETER :: h2ocoeffs(5, n_ace) = RESHAPE([ &
                                                       1.8275E-06, 1.9385E-06, 2.7918E-07, 1.1339E-06, 5.6496E-07, &
                                                       1.7482E-06, 1.9556E-06, 1.5120E-07, 1.2145E-06, 6.0461E-07, &
                                                       1.9014E-06, 1.8690E-06, 4.0900E-07, 1.0202E-06, 5.0890E-07, &
                                                       1.9885E-06, 1.7454E-06, 6.5235E-07, 9.1125E-07, 4.5717E-07, &
                                                       2.2510E-06, 1.9395E-06, -5.1958E-07, 6.0419E-07, 2.8430E-07, &
                                                       2.1314E-06, 1.6389E-06, -1.8774E-07, 9.6654E-07, 4.7675E-07, &
                                                       1.9797E-06, 1.5062E-06, 9.0113E-08, 8.8843E-07, 4.4129E-07, &
                                                       1.9711E-06, 1.2862E-06, 5.1189E-07, 5.0692E-07, 2.5295E-07, &
                                                       1.8712E-06, 1.1365E-06, 5.9111E-07, 3.5170E-07, 1.7584E-07, &
                                                       1.7060E-06, 1.0002E-06, 6.0783E-07, 3.2339E-07, 1.6288E-07, &
                                                       1.6449E-06, 6.3478E-07, 5.4028E-07, 3.3542E-07, 1.7185E-07, &
                                                       1.4023E-06, 5.5732E-07, 3.3814E-07, 3.1593E-07, 1.6028E-07, &
                                                       1.1985E-06, 4.3654E-07, 8.5842E-08, 4.5344E-07, 2.2897E-07, &
                                                       1.1847E-06, 2.5249E-07, -3.7605E-08, 5.8919E-07, 2.9888E-07, &
                                                       1.1189E-06, 2.5122E-07, 6.7708E-09, 5.6701E-07, 2.8805E-07, &
                                                       1.0350E-06, 6.5658E-08, -5.3691E-08, 5.7691E-07, 2.9424E-07, &
                                                       9.6944E-07, -1.7360E-08, -5.1035E-08, 5.0968E-07, 2.6065E-07, &
                                                       1.0475E-06, -2.0939E-08, -5.5097E-08, 5.8445E-07, 2.9894E-07, &
                                                       1.0109E-06, 1.1156E-07, 2.8096E-09, 6.2913E-07, 3.2120E-07, &
                                                       1.0611E-06, 3.2325E-08, 1.6240E-08, 5.7072E-07, 2.9222E-07, &
                                                       1.1450E-06, -8.9371E-09, 2.1557E-08, 6.0580E-07, 3.1066E-07, &
                                                       1.1811E-06, -1.5969E-07, 2.0227E-08, 5.8695E-07, 3.0247E-07, &
                                                       1.3049E-06, -3.4595E-07, 1.0663E-07, 6.0265E-07, 3.1336E-07, &
                                                       1.4610E-06, -3.9396E-07, -9.9152E-08, 4.6947E-07, 2.4321E-07, &
                                                       1.6839E-06, -5.2825E-07, -3.5202E-07, 3.3374E-07, 1.7206E-07, &
                                                       1.7722E-06, -7.4279E-07, -5.9396E-07, 1.2294E-07, 6.3375E-08, &
                                                       1.7731E-06, -9.1753E-07, -6.8094E-07, 2.3005E-07, 1.1895E-07, &
                                                       1.9041E-06, -1.2040E-06, -7.6971E-07, 2.6305E-07, 1.3765E-07, &
                                                       2.0364E-06, -1.3966E-06, -6.2757E-07, 2.4001E-07, 1.2940E-07, &
                                                       1.9232E-06, -1.6701E-06, -1.2561E-07, 5.7363E-07, 3.0884E-07, &
                                                       1.8864E-06, -1.8311E-06, 1.2302E-07, 7.5913E-07, 4.0835E-07, &
                                                       2.0659E-06, -1.8734E-06, 5.7529E-07, 3.9619E-07, 2.2813E-07, &
                                                       1.7182E-06, -1.5548E-06, -7.4811E-07, 8.4559E-07, 4.3975E-07, &
                                                       1.7492E-06, -1.6212E-06, -8.4923E-07, 7.3101E-07, 3.8054E-07, &
                                                       1.5048E-06, -1.7416E-06, -4.1964E-07, 1.0097E-06, 5.2947E-07, &
                                                       1.4903E-06, -1.7519E-06, -3.9055E-07, 1.0293E-06, 5.3995E-07], &
                                                       [5, n_ace])

!ACE-FTS 3-monthly O3 dataset.
!Fourier transform coefficients to 2nd order to make climatology.
!Olaf Morgenstern, Sept 2018
      REAL, PARAMETER :: o3coeffs(5, n_ace) = RESHAPE([ &
                                                      5.2093E-07, -3.0501E-07, 2.7602E-07, -1.3353E-07, -6.2171E-08, &
                                                      4.4081E-07, -2.6863E-07, 1.3607E-07, -4.4413E-08, -1.8514E-08, &
                                                      4.0440E-07, -2.3802E-07, 3.9693E-08, -6.1140E-09, -3.2167E-10, &
                                                      3.7188E-07, -1.6954E-07, -1.9669E-08, -2.4127E-08, -1.0916E-08, &
                                                      4.0728E-07, -2.4002E-07, 1.6186E-07, -1.6336E-08, -4.1134E-09, &
                                                      3.7676E-07, -1.4354E-07, 1.4899E-07, -6.8827E-08, -3.2103E-08, &
                                                      3.6341E-07, -1.0716E-07, 6.6159E-08, -2.0992E-08, -8.9248E-09, &
                                                      3.6908E-07, -1.1438E-07, -3.2458E-08, 3.7764E-08, 2.0094E-08, &
                                                      3.5972E-07, -8.2912E-08, -4.9043E-08, 4.0317E-08, 2.0898E-08, &
                                                      4.0258E-07, -1.1774E-07, -7.2229E-08, 4.3080E-08, 2.2386E-08, &
                                                      3.5882E-07, -5.7688E-08, -7.8390E-08, -1.1117E-08, -6.0398E-09, &
                                                      4.3687E-07, -1.0382E-07, -7.3913E-08, -2.2041E-08, -1.1128E-08, &
                                                      4.7485E-07, -6.6653E-08, -1.2075E-08, -1.2933E-07, -6.5730E-08, &
                                                      4.2807E-07, 1.2086E-09, 1.3166E-08, -1.3935E-07, -7.1241E-08, &
                                                      4.3500E-07, -8.1355E-08, -8.3265E-09, -8.5678E-08, -4.3182E-08, &
                                                      4.7448E-07, -7.4928E-08, 5.6054E-09, -1.3205E-07, -6.6838E-08, &
                                                      4.5172E-07, -4.6306E-08, 3.1445E-08, -1.4693E-08, -6.7030E-09, &
                                                      4.2373E-07, -4.6250E-08, 4.5800E-08, -1.0232E-08, -4.2510E-09, &
                                                      4.4011E-07, -3.0652E-08, 3.8499E-08, -7.2965E-08, -3.6624E-08, &
                                                      3.5624E-07, -6.6562E-08, 1.1966E-07, -4.2304E-08, -1.9618E-08, &
                                                      4.3206E-07, -4.7455E-08, 4.9168E-09, -1.3392E-07, -6.8075E-08, &
                                                      3.6728E-07, -3.5805E-08, 8.4126E-08, -7.9432E-08, -3.9354E-08, &
                                                      3.5121E-07, -4.0872E-08, 5.4060E-08, -1.2650E-07, -6.3764E-08, &
                                                      3.6722E-07, -4.9058E-08, 4.7771E-08, -8.2817E-08, -4.1381E-08, &
                                                      3.2403E-07, 4.2926E-10, 6.7643E-08, -4.0692E-08, -2.0060E-08, &
                                                      3.4871E-07, 5.0233E-08, 1.1321E-07, 5.8957E-08, 3.1025E-08, &
                                                      3.7537E-07, 6.8456E-08, 7.4805E-08, 6.4409E-08, 3.3190E-08, &
                                                      3.5921E-07, 7.5955E-08, 6.9903E-08, 8.6319E-08, 4.4282E-08, &
                                                      3.6728E-07, 9.1783E-08, 4.0672E-08, 1.2052E-07, 6.1304E-08, &
                                                      4.1543E-07, 1.2018E-07, -3.4157E-08, 8.7355E-08, 4.3163E-08, &
                                                      4.3857E-07, 1.3323E-07, -1.1552E-07, 3.1365E-08, 1.3406E-08, &
                                                      4.3702E-07, 1.8558E-07, -1.0531E-07, 6.7390E-08, 3.1462E-08, &
                                                      4.2658E-07, 1.3384E-07, 5.7944E-08, 5.5368E-08, 2.7717E-08, &
                                                      4.4652E-07, 1.7079E-07, 6.9623E-08, 7.5859E-08, 3.7985E-08, &
                                                      5.3677E-07, 2.5131E-07, 7.7937E-08, 7.1042E-08, 3.4821E-08, &
                                                      5.6287E-07, 2.6383E-07, 3.0677E-08, 4.1262E-08, 1.8892E-08], &
                                                      [5, n_ace])

      CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'UKCA_TOP_BOUNDARY'
      REAL(KIND=jprb)               :: zhook_handle
      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! make sure we initialize the array either if the years no longer match
! i.e. the job has crossed into a new year, of at the start of each leg
      IF ((firstcall) .OR. (i_year /= year_filled)) THEN
         ! determine the positions of the NO and CO tracers; this may need to be
         ! expanded.
         DO i = 1, jpctr
            IF (advt(i) == 'CO        ') i_co = i
            IF (advt(i) == 'NO        ') i_no = i
            IF (advt(i) == 'H2O       ') i_h2o = i
            IF (advt(i) == 'H2        ') i_h2 = i
            IF (advt(i) == 'H         ') i_h = i
            IF (advt(i) == 'OH        ') i_oh = i
            IF (advt(i) == 'O3        ') i_o3 = i
         END DO

         DO i = 1, n_ace
            acelat(i) = -90.0 + dellat_ace*(REAL(i) - 0.5)
         END DO

         ! Reconstitute field
         yearlen = days_in_year(i_year)
         ! allocate climatology of CO, NO, H2O, O3
         ALLOCATE (clim(n_ace, yearlen, n_top))

         ! expand climatology from Fourier coefficients; makes annually periodic field
         DO i = 1, n_ace
            DO j = 1, yearlen
               clim(i, j, j_co) = cocoeffs(1, i) + &
                                  cocoeffs(2, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  cocoeffs(3, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  cocoeffs(4, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi) + &
                                  cocoeffs(5, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi)
               clim(i, j, j_no) = nocoeffs(1, i) + &
                                  nocoeffs(2, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  nocoeffs(3, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  nocoeffs(4, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi) + &
                                  nocoeffs(5, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi)
               clim(i, j, j_h2o) = h2ocoeffs(1, i) + &
                                   h2ocoeffs(2, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                   h2ocoeffs(3, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                   h2ocoeffs(4, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi) + &
                                   h2ocoeffs(5, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi)
               clim(i, j, j_o3) = o3coeffs(1, i) + &
                                  o3coeffs(2, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  o3coeffs(3, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*2.0*pi) + &
                                  o3coeffs(4, i)*COS((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi) + &
                                  o3coeffs(5, i)*SIN((REAL(j) - 1.0)/REAL(yearlen)*4.0*pi)
            END DO
         END DO

         ! prevent negatives
         WHERE (clim < 10.0*peps) clim = 10.0*peps

         ! interpolate to UM grid
         IF (ALLOCATED(clim_interp)) DEALLOCATE (clim_interp)
         ALLOCATE (clim_interp(rows, yearlen, n_top))
         DO i = 1, rows
            frac = 0.0
            idx = FLOOR((latitude(i) - acelat(1))/dellat_ace) + 1
            IF ((idx >= 1) .AND. (idx < n_ace)) THEN
               frac = (latitude(i) - acelat(idx))/dellat_ace
               clim_interp(i, :, :) = (1.0 - frac)*clim(idx, :, :) + frac*clim(idx + 1, :, :)
            ELSE IF (idx < 1) THEN
               clim_interp(i, :, :) = clim(1, :, :)
            ELSE
               clim_interp(i, :, :) = clim(n_ace, :, :)
            END IF
         END DO

         ! prevent negatives
         WHERE (clim_interp < 10.0*peps) clim_interp = 10.0*peps

         ! rescale to make MMR
         clim_interp(:, :, j_no) = clim_interp(:, :, j_no)*c_no
         clim_interp(:, :, j_co) = clim_interp(:, :, j_co)*c_co
         clim_interp(:, :, j_h2o) = clim_interp(:, :, j_h2o)*c_h2o
         clim_interp(:, :, j_o3) = clim_interp(:, :, j_o3)*c_o3
         firstcall = .FALSE.
         year_filled = i_year
         DEALLOCATE (clim)
      END IF

! impose upper-boundary mixing ratio for NO, CO, O3, and H2O.
      DO i = 1, row_length
         tracer(i, :, model_levels, i_no) = clim_interp(:, i_day_number, j_no)
         tracer(i, :, model_levels, i_co) = clim_interp(:, i_day_number, j_co)
         tracer(i, :, model_levels, i_o3) = clim_interp(:, i_day_number, j_o3)
      END DO

! Take note of change of water vapour; reduce H2, H and OH accordingly to
! preserve hydrogen
      IF (ukca_config%i_ukca_topboundary == i_top_BC_H2O) THEN
         ALLOCATE (change_h2o(rows))
         ALLOCATE (toth(rows))
         DO i = 1, row_length
            change_h2o = (clim_interp(:, i_day_number, j_h2o) - &
                          tracer(i, :, model_levels, i_h2o))/c_h2o
            tracer(i, :, model_levels, i_h2o) = clim_interp(:, i_day_number, j_h2o)
            ! calculate total hydrogen for OH, H, H2 (the other significant hydrogen
            ! compounds at 85 km).
            toth = tracer(i, :, model_levels, i_oh)*0.5/c_oh &
                   + tracer(i, :, model_levels, i_h)*0.5/c_h &
                   + tracer(i, :, model_levels, i_h2)/c_h2
            ! calculate rescaling factor for these compounds
            toth = MAX(1.0 - change_h2o/MAX(toth, 1.0E-20), 0.0)
            tracer(i, :, model_levels, i_h2) = tracer(i, :, model_levels, i_h2)*toth
            tracer(i, :, model_levels, i_h) = tracer(i, :, model_levels, i_h)*toth
            tracer(i, :, model_levels, i_oh) = tracer(i, :, model_levels, i_oh)*toth
         END DO
         DEALLOCATE (change_h2o)
         DEALLOCATE (toth)
      END IF

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

   END SUBROUTINE ukca_top_boundary

END MODULE ukca_topboundary_mod
