! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Subroutine to calculate impaction scavenging of aerosols
!    by falling raindrops. This is a two-moment scavenging scheme
!    for liquid raindroplets only, and currently only used
!    for insoluble dust
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ######################################################################
!
! .. Subroutine Interface:
MODULE ukca_impc_scav_dust_mod

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim

   IMPLICIT NONE

   PRIVATE

! .. Variables for impaction scavenging
   INTEGER, PARAMETER  :: ncol = 22 ! # of aero diameter integration pts
   INTEGER, PARAMETER  :: nrow = 22 ! # of rainfall rate integration pts
   INTEGER, PARAMETER  :: nlev = 5  ! # of aero stdev integration pts

! .. Parameters for st. dev. integration points
! .. Formulae: (i-1)*dt_stdev + stdev_low for i=1,..,nlev
   REAL, PARAMETER :: stdev_low = 1.2
   REAL, PARAMETER :: dt_stdev = 0.2

! .. Parameters for median diameter integration points
! .. Formulae: 2 * 10**((i-1)*dt_ldp + ldplow) for i=1,..,ncol
   REAL, PARAMETER :: ldplow = -9.0
   REAL, PARAMETER :: dt_ldp = 0.2

! .. Parameters to define rainrate integration points
! .. Formulate: 10**((i-1)*dt_lrf + lrflow) for i=1,..,nrow
   REAL, PARAMETER :: lrflow = -1.0
   REAL, PARAMETER :: dt_lrf = 1.0/7.0

! .. Scavenging  coefficients (s-1) for aerosol number and mass
! .. as a function of median diameter, rainrate and modal stdev
! .. Populated in ukca_impc_scav_dust_init() below
   REAL, ALLOCATABLE, SAVE :: scav_coeff_num(:, :, :)
   REAL, ALLOCATABLE, SAVE :: scav_coeff_mass(:, :, :)

! Routines available from this module
   PUBLIC :: ukca_impc_scav_dust, ukca_impc_scav_dust_init, &
             ukca_impc_scav_dust_dealloc

   CHARACTER(LEN=*), PARAMETER :: ModuleName = 'UKCA_IMPC_SCAV_DUST_MOD'

CONTAINS

! ----------------------------------------------------------------------
   SUBROUTINE ukca_impc_scav_dust(nbox, nbudaer, nd, md, mdt, crain, drain, &
                                  wetdp, dtc, bud_aer_mas, iextra_checks)
! ----------------------------------------------------------------------
!
! Purpose:
! -------
!     Subroutine to calculate impaction scavenging of aerosols
!     by falling raindrops.
!
!     The method differs to ukca_impc_scav and depends on the values
!     in the scav_coeff_num and scav_coeff_mass arrays. Scavenging
!     coefficients are either derived from Slinn (1983) formulae
!     for collision efficiency, integrated over the rain droplet
!     size distribution (Abel and Boutle, 2012) and using droplet
!     fall speeds from Beard (1976) or provided directly from a
!     parameterisation in Laakso et al (2003)
!
!     The rainfall rate and the geometric median diameters integration
!     points are dispersed logarithmically while the standard deviation
!     integration points are (1.2, 1.4, ..., 2). Interpolation of the
!     scavenging coefficient arrays is log10-log10 with diameter (i.e.
!     that both diameter and scav coefficient are interpolated in log10
!     space) and linear with rain rate, and using a nearest-neighbour
!     approach for standard deviation
!
!     Parameters
!     ----------
!
!     Inputs
!     ------
!     NBOX        : Number of grid boxes
!     NBUDAER     : Number of aerosol budget fields
!     ND          : Aerosol ptcl number density (cm^-3)
!     MD          : Avg cpt mass of aerosol ptcl (particle^-1)
!     MDT         : Total median aerosol mass (molecules per ptcl)
!     CRAIN       : Convective rain rate array (kgm^-2s^-1)
!     DRAIN       : Dynamic rain rate array (kgm^-2s^-1)
!     WETDP       : Wet diameter corresponding to DRYDP (m)
!     DTC         : Time step of process (s)
!
!     Outputs
!     -------
!     ND          : new aerosol number conc (cm^-3)
!     MD          : new cpt mass of aerosol ptcl (particle^-1)
!     MDT         : Total median aerosol mass (molecules per ptcl)
!     BUD_AER_MAS : Updated aerosol budgets
!
!----------------------------------------------------------------------

      USE ukca_config_specification_mod, ONLY: glomap_variables

      USE ukca_mode_setup, ONLY: cp_du, nmodes, mode_acc_insol, mode_cor_insol, &
                                 mode_sup_insol, cp_mp

      USE ukca_setup_indices, ONLY: nmasimscduaccins, nmasimscducorins, &
                                    nmasimscdusupins, nmasimscmpaccins, &
                                    nmasimscmpcorins, nmasimscmpsupins

      USE ukca_mode_check_artefacts_mod, ONLY: ukca_mode_check_mdt
      USE ukca_types_mod, ONLY: logical_32

      IMPLICIT NONE

! .. Subroutine interface
      INTEGER, INTENT(IN)  :: nbox
      INTEGER, INTENT(IN)  :: nbudaer
      INTEGER, INTENT(IN)  :: iextra_checks
      REAL, INTENT(IN)     :: wetdp(nbox, nmodes)
      REAL, INTENT(IN)     :: dtc
      REAL, INTENT(IN)     :: crain(nbox)
      REAL, INTENT(IN)     :: drain(nbox)
      REAL, INTENT(IN OUT) :: nd(nbox, nmodes)
      REAL, INTENT(IN OUT) :: md(nbox, nmodes, glomap_variables%ncp)
      REAL, INTENT(IN OUT) :: mdt(nbox, nmodes)
      REAL, INTENT(IN OUT) :: bud_aer_mas(nbox, 0:nbudaer)

! .. Local variables

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
      LOGICAL, POINTER :: component(:, :)
      REAL, POINTER :: mfrac_0(:, :)
      REAL, POINTER :: mmid(:)
      LOGICAL, POINTER :: mode(:)
      INTEGER, POINTER :: ncp
      REAL, POINTER :: num_eps(:)
      REAL, POINTER :: sigmag(:)

      REAL, PARAMETER :: fc = 0.3
      REAL, PARAMETER :: fd = 1.0
      REAL, PARAMETER :: secs_per_hr = 3600.0
      REAL    :: allfrac(2)
      REAL    :: allrain(2, nbox)
      REAL    :: totrain(nbox)
      REAL    :: scavn(2, nbox, nmodes)
      REAL    :: scavm(2, nbox, nmodes)
      REAL    :: deln, deln1, deln2
      REAL    :: dm1(glomap_variables%ncp)
      REAL    :: dm2(glomap_variables%ncp)
      REAL    :: dm(glomap_variables%ncp)
      REAL    :: ndnew
      REAL    :: ii, dplow, dpupp, jj, rflow, rfupp
      REAL    :: logdp, logdplow, logdpupp, fac1, fac2
      REAL    :: sc_num_dp_rf1, sc_num_dp_rf2
      REAL    :: sc_mass_dp_rf1, sc_mass_dp_rf2
      LOGICAL :: l_interp_dp, l_interp_RF
      LOGICAL(KIND=logical_32) :: mask(nbox)
      INTEGER :: ilow, iupp, jlow, jupp, k
      INTEGER :: imode, icp, jl, iprecip

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_IMPC_SCAV_DUST'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Caution - pointers to TYPE glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
      component => glomap_variables%component
      mfrac_0 => glomap_variables%mfrac_0
      mmid => glomap_variables%mmid
      mode => glomap_variables%mode
      ncp => glomap_variables%ncp
      num_eps => glomap_variables%num_eps
      sigmag => glomap_variables%sigmag

! .. Combine the convective and dynamic rain in an array
! .. to loop over and convert rain from kgm-2s-1 to mm/hr
      allfrac(1) = fc
      allfrac(2) = fd
      allrain(:, :) = 0.0
      totrain(:) = 0.0
      DO jl = 1, nbox
         allrain(1, jl) = crain(jl)*secs_per_hr/allfrac(1)
         allrain(2, jl) = drain(jl)*secs_per_hr/allfrac(2)
         totrain(jl) = (crain(jl)*secs_per_hr) + (drain(jl)*secs_per_hr)
      END DO

! .. Initialise 3D scavenging arrays
      scavn(:, :, :) = 0.0
      scavm(:, :, :) = 0.0

! .. Loop over mode, rain type, and gridcell to populate the
! .. 3D real-time scavenging arrays
      DO imode = mode_acc_insol, mode_sup_insol
         IF (mode(imode)) THEN

            ! .. Nearest neighbour interpolation is used for standard
            ! .. deviation in scavenging arrays so only need to do
            ! .. this once for each mode
            k = MAX(NINT((sigmag(imode) - stdev_low + dt_stdev)/dt_stdev), 1)
            k = MIN(k, nlev)

            DO jl = 1, nbox
               DO iprecip = 1, 2
                  IF (allrain(iprecip, jl) > 0.0) THEN

                     ! .. Find the i index (diameter) for interpolation of
                     ! .. scavenging coefficient arrays
                     l_interp_dp = .FALSE.
                     ii = 1 + (LOG10(wetdp(jl, imode)/2.0) - ldplow)/dt_ldp
                     IF (ii <= 1) THEN
                        ilow = 1
                        iupp = 1
                        dplow = 2.0*(10.0**ldplow)
                        dpupp = dplow
                     ELSE IF (ii >= ncol) THEN
                        ilow = ncol
                        iupp = ncol
                        dplow = 2.0*(10.0**(ldplow + (ilow - 1)*dt_ldp))
                        dpupp = dplow
                     ELSE
                        ilow = INT(ii)
                        iupp = ilow + 1
                        dplow = 2.0*(10.0**(ldplow + (ilow - 1)*dt_ldp))
                        dpupp = 2.0*(10.0**(ldplow + (iupp - 1)*dt_ldp))
                        l_interp_dp = .TRUE.
                     END IF

                     ! .. Find the j index (rainfall) for interpolation of
                     ! .. scavenging coefficient arrays
                     l_interp_RF = .FALSE.
                     jj = 1 + (LOG10(allrain(iprecip, jl)) - lrflow)/dt_lrf
                     IF (jj <= 1) THEN
                        jlow = 1
                        jupp = 1
                        rflow = 10.0**lrflow
                        rfupp = rflow
                     ELSE IF (jj >= nrow) THEN
                        jlow = nrow
                        jupp = nrow
                        rflow = 10.0**(lrflow + (jlow - 1)*dt_lrf)
                        rfupp = rflow
                     ELSE
                        jlow = INT(jj)
                        jupp = jlow + 1
                        rflow = 10.0**(lrflow + (jlow - 1)*dt_lrf)
                        rfupp = 10.0**(lrflow + (jupp - 1)*dt_lrf)
                        l_interp_RF = .TRUE.
                     END IF

                     ! .. Now interpolate in i index (LOG10-LOG10)
                     IF (l_interp_dp) THEN
                        logdp = LOG10(wetdp(jl, imode))
                        logdplow = LOG10(dplow)
                        logdpupp = LOG10(dpupp)
                        fac1 = (logdpupp - logdp)/(logdpupp - logdplow)
                        fac2 = (logdp - logdplow)/(logdpupp - logdplow)
                        sc_num_dp_rf1 = &
                           10.0**(fac1*LOG10(scav_coeff_num(ilow, jlow, k)) + &
                                  fac2*LOG10(scav_coeff_num(iupp, jlow, k)))
                        sc_num_dp_rf2 = &
                           10.0**(fac1*LOG10(scav_coeff_num(ilow, jupp, k)) + &
                                  fac2*LOG10(scav_coeff_num(iupp, jupp, k)))
                        sc_mass_dp_rf1 = &
                           10.0**(fac1*LOG10(scav_coeff_mass(ilow, jlow, k)) + &
                                  fac2*LOG10(scav_coeff_mass(iupp, jlow, k)))
                        sc_mass_dp_rf2 = &
                           10.0**(fac1*LOG10(scav_coeff_mass(ilow, jupp, k)) + &
                                  fac2*LOG10(scav_coeff_mass(iupp, jupp, k)))
                     ELSE
                        sc_num_dp_rf1 = scav_coeff_num(ilow, jlow, k)
                        sc_num_dp_rf2 = scav_coeff_num(ilow, jupp, k)
                        sc_mass_dp_rf1 = scav_coeff_mass(ilow, jlow, k)
                        sc_mass_dp_rf2 = scav_coeff_mass(ilow, jupp, k)
                     END IF

                     ! .. Lastly interpolate in the j index (LIN-LIN)
                     IF (l_interp_RF) THEN
                        fac1 = (rfupp - allrain(iprecip, jl))/(rfupp - rflow)
                        fac2 = (allrain(iprecip, jl) - rflow)/(rfupp - rflow)
                        scavn(iprecip, jl, imode) = fac1*sc_num_dp_rf1 + &
                                                    fac2*sc_num_dp_rf2
                        scavm(iprecip, jl, imode) = fac1*sc_mass_dp_rf1 + &
                                                    fac2*sc_mass_dp_rf2
                     ELSE
                        scavn(iprecip, jl, imode) = sc_num_dp_rf1
                        scavm(iprecip, jl, imode) = sc_mass_dp_rf1
                     END IF

                  END IF ! .. IF (allrain(iprecip,jl) > 0.0)
               END DO ! .. iprecip
            END DO ! .. nbox
         END IF ! .. IF (mode(imode))
      END DO ! .. imode

! .. Apply derived scavenging coefficients to input mass/number
      DO jl = 1, nbox
         IF (totrain(jl) > 0.0) THEN
            DO imode = mode_acc_insol, mode_sup_insol
               IF (mode(imode)) THEN

                  ! .. Only do anything if the initial number is greater than
                  ! .. a threshold value (num_eps)
                  IF (nd(jl, imode) > num_eps(imode)) THEN

                     ! .. Nullify mass and number changes
                     DO icp = 1, ncp
                        dm1(icp) = 0.0
                        dm2(icp) = 0.0
                     END DO
                     deln1 = 0.0
                     deln2 = 0.0

                     ! .. Convective rain
                     IF (crain(jl) > 0.0) THEN
                        deln1 = allfrac(1)*nd(jl, imode)*(1.0 - EXP(-scavn(1, jl, imode)*dtc))
                        DO icp = 1, ncp
                           IF (component(imode, icp)) THEN
                              dm1(icp) = allfrac(1)*nd(jl, imode)*md(jl, imode, icp)* &
                                         (1.0 - EXP(-scavm(1, jl, imode)*dtc))
                           END IF
                        END DO
                     END IF

                     ! .. Dynamical rain
                     IF (drain(jl) > 0.0) THEN
                        deln2 = allfrac(2)*nd(jl, imode)*(1.0 - EXP(-scavn(2, jl, imode)*dtc))
                        DO icp = 1, ncp
                           IF (component(imode, icp)) THEN
                              dm2(icp) = allfrac(2)*nd(jl, imode)*md(jl, imode, icp)* &
                                         (1.0 - EXP(-scavm(2, jl, imode)*dtc))
                           END IF
                        END DO
                     END IF

                     ! .. Sum mass and number changes from dynamic and convective rain
                     deln = MIN(nd(jl, imode), deln1 + deln2)
                     DO icp = 1, ncp
                        IF (component(imode, icp)) THEN
                           dm(icp) = MIN(nd(jl, imode)*md(jl, imode, icp), dm1(icp) + dm2(icp))
                        END IF
                     END DO

                     ! .. Update number and mass concentrations
                     mdt(jl, imode) = 0.0
                     ndnew = nd(jl, imode) - deln
                     IF (ndnew > num_eps(imode)) THEN
                        DO icp = 1, ncp
                           IF (component(imode, icp)) THEN
                              md(jl, imode, icp) = (nd(jl, imode)*md(jl, imode, icp) - dm(icp))/ndnew
                              mdt(jl, imode) = mdt(jl, imode) + md(jl, imode, icp)
                           END IF
                        END DO
                        nd(jl, imode) = ndnew
                     ELSE
                        deln = nd(jl, imode)
                        DO icp = 1, ncp
                           IF (component(imode, icp)) THEN
                              dm(icp) = nd(jl, imode)*md(jl, imode, icp)
                              md(jl, imode, icp) = mmid(imode)*mfrac_0(imode, icp)
                           END IF
                        END DO
                        nd(jl, imode) = 0.0
                        mdt(jl, imode) = mmid(imode)
                     END IF

                     ! .. Store cpt imp scav mass fluxes for budget calculations
                     DO icp = 1, ncp
                        IF (component(imode, icp)) THEN
                           IF (icp == cp_du) THEN
                              IF ((imode == mode_acc_insol) .AND. (nmasimscduaccins > 0)) &
                                 bud_aer_mas(jl, nmasimscduaccins) = &
                                 bud_aer_mas(jl, nmasimscduaccins) + dm(icp)
                              IF ((imode == mode_cor_insol) .AND. (nmasimscducorins > 0)) &
                                 bud_aer_mas(jl, nmasimscducorins) = &
                                 bud_aer_mas(jl, nmasimscducorins) + dm(icp)
                              IF ((imode == mode_sup_insol) .AND. (nmasimscdusupins > 0)) &
                                 bud_aer_mas(jl, nmasimscdusupins) = &
                                 bud_aer_mas(jl, nmasimscdusupins) + dm(icp)
                           END IF
                           IF (icp == cp_mp) THEN
                              IF ((imode == mode_acc_insol) .AND. (nmasimscmpaccins > 0)) &
                                 bud_aer_mas(jl, nmasimscmpaccins) = &
                                 bud_aer_mas(jl, nmasimscmpaccins) + dm(icp)
                              IF ((imode == mode_cor_insol) .AND. (nmasimscmpcorins > 0)) &
                                 bud_aer_mas(jl, nmasimscmpcorins) = &
                                 bud_aer_mas(jl, nmasimscmpcorins) + dm(icp)
                              IF ((imode == mode_sup_insol) .AND. (nmasimscmpsupins > 0)) &
                                 bud_aer_mas(jl, nmasimscmpsupins) = &
                                 bud_aer_mas(jl, nmasimscmpsupins) + dm(icp)
                           END IF
                        END IF ! .. if component present
                     END DO ! .. icp

                  END IF ! .. IF (nd(jl,imode) > num_eps(imode))
               END IF ! .. IF (mode(imode))
            END DO ! .. imode
         END IF !.. IF (totrain(jl) > 0.0)
      END DO ! .. jl

! Apply extra checks on MDT being out of range after impaction scavenging
      IF (iextra_checks > 1) THEN !do this only when ie_c = 2
         DO imode = mode_acc_insol, mode_sup_insol
            IF (mode(imode)) THEN
               mask(:) = .TRUE. ! dummy mask
               CALL ukca_mode_check_mdt(nbox, imode, mdt, md, nd, mask)
            END IF
         END DO
      END IF

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_impc_scav_dust

! ----------------------------------------------------------------------
   SUBROUTINE ukca_impc_scav_dust_init(verbose)
! ----------------------------------------------------------------------
! Purpose:
! -------
! Initialise scavenging coefficient lookup tables
! ----------------------------------------------------------------------

      USE umPrintMgr, ONLY: umPrint, umMessage

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: verbose ! flag to indicate level of verbosity

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_IMPC_SCAV_DUST_INIT'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! New dust impaction scavenging scheme - this module
      IF (.NOT. ALLOCATED(scav_coeff_num)) &
         ALLOCATE (scav_coeff_num(1:ncol, 1:nrow, 1:nlev))
      IF (.NOT. ALLOCATED(scav_coeff_mass)) &
         ALLOCATE (scav_coeff_mass(1:ncol, 1:nrow, 1:nlev))

      IF (verbose >= 2) THEN
         WRITE (umMessage, '(A50,2I5)') 'New dust impaction scavenging scheme on '// &
            'in this setup'
         CALL umPrint(umMessage, src=RoutineName)
      END IF

! Populate the new scavenging arrays
      scav_coeff_num(1:ncol, 1, 1) = [ &
                                     6.1037148E-05, 3.4363694E-05, 1.9950263E-05, 1.1905033E-05, 7.2960632E-06, &
                                     4.6040514E-06, 3.0113628E-06, 2.0628355E-06, 1.4976336E-06, 1.1632750E-06, &
                                     9.7163390E-07, 8.7728381E-07, 8.6524578E-07, 9.4244138E-07, 1.1398267E-06, &
                                     1.6519928E-06, 1.2224743E-05, 4.3962760E-05, 5.5330641E-05, 5.4394218E-05, &
                                     6.0671539E-05, 9.0363565E-05]
      scav_coeff_num(1:ncol, 1, 2) = [ &
                                     6.5706309E-05, 3.6689561E-05, 2.1157495E-05, 1.2556026E-05, 7.6590291E-06, &
                                     4.8121111E-06, 3.1332367E-06, 2.1353433E-06, 1.5412965E-06, 1.1902741E-06, &
                                     9.9004062E-07, 8.9285139E-07, 8.8212565E-07, 9.6559994E-07, 1.2355951E-06, &
                                     3.4073185E-06, 1.6383684E-05, 4.0091787E-05, 5.3066348E-05, 5.5751309E-05, &
                                     6.5116983E-05, 1.0058589E-04]
      scav_coeff_num(1:ncol, 1, 3) = [ &
                                     7.2835969E-05, 4.0195060E-05, 2.2954348E-05, 1.3514131E-05, 8.1881546E-06, &
                                     5.1130763E-06, 3.3084738E-06, 2.2391734E-06, 1.6037626E-06, 1.2290530E-06, &
                                     1.0165489E-06, 9.1507999E-07, 9.0711056E-07, 1.0318727E-06, 1.7965781E-06, &
                                     6.0411697E-06, 1.8899684E-05, 3.7585026E-05, 5.0641188E-05, 5.7479621E-05, &
                                     7.1931053E-05, 1.1693697E-04]
      scav_coeff_num(1:ncol, 1, 4) = [ &
                                     8.2525145E-05, 4.4882987E-05, 2.5319472E-05, 1.4756997E-05, 8.8659539E-06, &
                                     5.4946394E-06, 3.5288539E-06, 2.3690317E-06, 1.6817423E-06, 1.2775862E-06, &
                                     1.0498017E-06, 9.4370470E-07, 9.5525410E-07, 1.2631306E-06, 2.8580057E-06, &
                                     8.4259391E-06, 2.0457544E-05, 3.6127692E-05, 4.9081412E-05, 5.9891848E-05, &
                                     8.1328636E-05, 1.4037703E-04]
      scav_coeff_num(1:ncol, 1, 5) = [ &
                                     9.5177508E-05, 5.0897106E-05, 2.8299513E-05, 1.6296551E-05, 9.6929965E-06, &
                                     5.9543943E-06, 3.7917755E-06, 2.5228860E-06, 1.7738382E-06, 1.3349742E-06, &
                                     1.0898021E-06, 9.8577057E-07, 1.0692644E-06, 1.7122465E-06, 4.1157658E-06, &
                                     1.0357743E-05, 2.1534589E-05, 3.5427627E-05, 4.8691800E-05, 6.3622664E-05, &
                                     9.4055582E-05, 1.7273542E-04]
      scav_coeff_num(1:ncol, 2, 1) = [ &
                                     6.3630249E-05, 3.5905266E-05, 2.0864780E-05, 1.2451834E-05, 7.6307078E-06, &
                                     4.8191336E-06, 3.1615601E-06, 2.1794014E-06, 1.5964712E-06, 1.2498167E-06, &
                                     1.0454355E-06, 9.3824472E-07, 9.1931965E-07, 1.0024244E-06, 1.2237700E-06, &
                                     1.8837638E-06, 1.7546948E-05, 5.7798820E-05, 6.9506943E-05, 6.7048700E-05, &
                                     7.3235012E-05, 1.0653566E-04]
      scav_coeff_num(1:ncol, 2, 2) = [ &
                                     6.8425657E-05, 3.8311325E-05, 2.2121036E-05, 1.3132244E-05, 8.0111851E-06, &
                                     5.0375335E-06, 3.2893441E-06, 2.2549185E-06, 1.6411685E-06, 1.2767445E-06, &
                                     1.0638332E-06, 9.5493724E-07, 9.3861199E-07, 1.0292958E-06, 1.3619858E-06, &
                                     4.4893250E-06, 2.2390521E-05, 5.2423241E-05, 6.6889678E-05, 6.8672629E-05, &
                                     7.8383078E-05, 1.1833013E-04]
      scav_coeff_num(1:ncol, 2, 3) = [ &
                                     7.5728623E-05, 4.1929178E-05, 2.3987342E-05, 1.4132263E-05, 8.5653106E-06, &
                                     5.3532468E-06, 3.4730013E-06, 2.3630820E-06, 1.7052896E-06, 1.3157713E-06, &
                                     1.0906415E-06, 9.7873724E-07, 9.6731271E-07, 1.1188625E-06, 2.1728461E-06, &
                                     8.0682766E-06, 2.5227189E-05, 4.8777900E-05, 6.3873277E-05, 7.0709515E-05, &
                                     8.6274862E-05, 1.3719369E-04]
      scav_coeff_num(1:ncol, 2, 4) = [ &
                                     8.5620121E-05, 4.6752579E-05, 2.6437645E-05, 1.5426999E-05, 9.2741656E-06, &
                                     5.7531346E-06, 3.7038482E-06, 2.4984124E-06, 1.7855925E-06, 1.3650462E-06, &
                                     1.1246306E-06, 1.0098031E-06, 1.0283858E-06, 1.4431214E-06, 3.6365487E-06, &
                                     1.1170250E-05, 2.6928043E-05, 4.6565478E-05, 6.1758467E-05, 7.3460738E-05, &
                                     9.7141839E-05, 1.6423156E-04]
      scav_coeff_num(1:ncol, 2, 5) = [ &
                                     9.8488088E-05, 5.2918404E-05, 2.9515462E-05, 1.7026882E-05, 1.0137556E-05, &
                                     6.2343828E-06, 3.9790803E-06, 2.6588234E-06, 1.8807314E-06, 1.4237662E-06, &
                                     1.1661790E-06, 1.0587062E-06, 1.1815658E-06, 2.0598127E-06, 5.3174238E-06, &
                                     1.3617522E-05, 2.8071647E-05, 4.5370769E-05, 6.0995795E-05, 7.7665489E-05, &
                                     1.1182145E-04, 2.0154434E-04]
      scav_coeff_num(1:ncol, 3, 1) = [ &
                                     6.6764925E-05, 3.7768020E-05, 2.1972740E-05, 1.3118339E-05, 8.0435134E-06, &
                                     5.0905131E-06, 3.3582677E-06, 2.3395106E-06, 1.7383607E-06, 1.3775219E-06, &
                                     1.1545280E-06, 1.0251831E-06, 9.9051368E-07, 1.0748815E-06, 1.3205944E-06, &
                                     2.1970309E-06, 2.4881232E-05, 7.5190835E-05, 8.7106001E-05, 8.2696223E-05, &
                                     8.8475609E-05, 1.2553365E-04]
      scav_coeff_num(1:ncol, 3, 2) = [ &
                                     7.1717208E-05, 4.0271805E-05, 2.3288238E-05, 1.3834254E-05, 8.4452229E-06, &
                                     5.3215353E-06, 3.4932499E-06, 2.4185050E-06, 1.7838321E-06, 1.4036331E-06, &
                                     1.1722734E-06, 1.0430865E-06, 1.0130085E-06, 1.1066022E-06, 1.5214396E-06, &
                                     5.9807808E-06, 3.0316801E-05, 6.8018367E-05, 8.4075898E-05, 8.4608445E-05, &
                                     9.4414322E-05, 1.3908811E-04]
      scav_coeff_num(1:ncol, 3, 3) = [ &
                                     7.9237922E-05, 4.4027369E-05, 2.5238717E-05, 1.4884919E-05, 9.0296396E-06, &
                                     5.6552044E-06, 3.6871153E-06, 2.5316617E-06, 1.8493319E-06, 1.4420764E-06, &
                                     1.1986962E-06, 1.0685709E-06, 1.0466121E-06, 1.2291029E-06, 2.6768051E-06, &
                                     1.0771180E-05, 3.3433724E-05, 6.2938711E-05, 8.0352298E-05, 8.6977394E-05, &
                                     1.0351285E-04, 1.6076000E-04]
      scav_coeff_num(1:ncol, 3, 4) = [ &
                                     8.9388215E-05, 4.9018315E-05, 2.7792733E-05, 1.6242450E-05, 9.7761106E-06, &
                                     6.0773440E-06, 3.9305950E-06, 2.6732936E-06, 1.9317611E-06, 1.4913466E-06, &
                                     1.2328316E-06, 1.1024455E-06, 1.1254178E-06, 1.6804280E-06, 4.6703803E-06, &
                                     1.4765215E-05, 3.5245300E-05, 5.9751670E-05, 7.7543680E-05, 9.0085815E-05, &
                                     1.1601922E-04, 1.9181368E-04]
      scav_coeff_num(1:ncol, 3, 5) = [ &
                                     1.0253986E-04, 5.5374294E-05, 3.0990421E-05, 1.7915601E-05, 1.0683545E-05, &
                                     6.5846386E-06, 4.2206325E-06, 2.8412692E-06, 2.0298892E-06, 1.5508223E-06, &
                                     1.2756222E-06, 1.1601750E-06, 1.3315458E-06, 2.5194053E-06, 6.8952972E-06, &
                                     1.7843860E-05, 3.6426548E-05, 5.7903430E-05, 7.6285200E-05, 9.4786815E-05, &
                                     1.3287082E-04, 2.3465061E-04]
      scav_coeff_num(1:ncol, 4, 1) = [ &
                                     7.0469649E-05, 3.9968562E-05, 2.3284417E-05, 1.3911934E-05, 8.5413885E-06, &
                                     5.4264945E-06, 3.6126266E-06, 2.5579084E-06, 1.9412913E-06, 1.5657118E-06, &
                                     1.3163211E-06, 1.1505808E-06, 1.0857031E-06, 1.1628835E-06, 1.4318934E-06, &
                                     2.6245843E-06, 3.4781982E-05, 9.6934846E-05, 1.0895754E-04, 1.0205357E-04, &
                                     1.0698456E-04, 1.4786945E-04]
      scav_coeff_num(1:ncol, 4, 2) = [ &
                                     7.5611234E-05, 4.2588656E-05, 2.4670026E-05, 1.4669901E-05, 8.9684048E-06, &
                                     5.6726922E-06, 3.7562388E-06, 2.6407795E-06, 1.9869166E-06, 1.5896531E-06, &
                                     1.3322056E-06, 1.1696531E-06, 1.1125274E-06, 1.2011142E-06, 1.7234466E-06, &
                                     8.0071948E-06, 4.0650586E-05, 8.7636738E-05, 1.0542731E-04, 1.0427307E-04, &
                                     1.1381060E-04, 1.6338662E-04]
      scav_coeff_num(1:ncol, 4, 3) = [ &
                                     8.3396747E-05, 4.6508803E-05, 2.6720336E-05, 1.5780583E-05, 9.5888859E-06, &
                                     6.0278768E-06, 3.9622694E-06, 2.7594745E-06, 2.0530482E-06, 1.6259471E-06, &
                                     1.3569096E-06, 1.1967726E-06, 1.1526806E-06, 1.3696665E-06, 3.3459956E-06, &
                                     1.4333410E-05, 4.3985798E-05, 8.0771823E-05, 1.0084737E-04, 1.0699049E-04, &
                                     1.2425512E-04, 1.8818449E-04]
      scav_coeff_num(1:ncol, 4, 4) = [ &
                                     9.3865615E-05, 5.1701287E-05, 2.9397774E-05, 1.7212622E-05, 1.0380101E-05, &
                                     6.4765864E-06, 4.2207127E-06, 2.9080839E-06, 2.1368933E-06, 1.6737241E-06, &
                                     1.3899735E-06, 1.2337282E-06, 1.2556044E-06, 1.9922748E-06, 6.0297663E-06, &
                                     1.9430968E-05, 4.5863510E-05, 7.6342969E-05, 9.7169834E-05, 1.1046360E-04, &
                                     1.3858109E-04, 2.2370008E-04]
      scav_coeff_num(1:ncol, 4, 5) = [ &
                                     1.0737297E-04, 5.8288164E-05, 3.2738823E-05, 1.8972887E-05, 1.1339904E-05, &
                                     7.0148886E-06, 4.5282008E-06, 3.0844549E-06, 2.2374451E-06, 1.7326860E-06, &
                                     1.4331763E-06, 1.3026657E-06, 1.5324700E-06, 3.1229744E-06, 8.9489268E-06, &
                                     2.3281278E-05, 4.7043931E-05, 7.3644645E-05, 9.5258454E-05, 1.1567048E-04, &
                                     1.5783441E-04, 2.7266614E-04]
      scav_coeff_num(1:ncol, 5, 1) = [ &
                                     7.4777212E-05, 4.2526200E-05, 2.4811988E-05, 1.4841745E-05, 9.1334172E-06, &
                                     5.8386359E-06, 3.9407403E-06, 2.8564586E-06, 2.2324790E-06, 1.8441613E-06, &
                                     1.5581269E-06, 1.3343411E-06, 1.2159424E-06, 1.2712274E-06, 1.5598280E-06, &
                                     3.2084270E-06, 4.7897526E-05, 1.2400555E-04, 1.3608372E-04, 1.2600302E-04, &
                                     1.2948627E-04, 1.7415961E-04]
      scav_coeff_num(1:ncol, 5, 2) = [ &
                                     8.0142562E-05, 4.5282389E-05, 2.6279349E-05, 1.5648879E-05, 9.5902893E-06, &
                                     6.1029521E-06, 4.0946420E-06, 2.9435275E-06, 2.2771123E-06, 1.8636369E-06, &
                                     1.5700150E-06, 1.3542698E-06, 1.2487074E-06, 1.3184204E-06, 1.9800422E-06, &
                                     1.0719179E-05, 5.3973254E-05, 1.1220209E-04, 1.3193118E-04, 1.2854522E-04, &
                                     1.3730564E-04, 1.9185783E-04]
      scav_coeff_num(1:ncol, 5, 3) = [ &
                                     8.8242873E-05, 4.9395731E-05, 2.8446256E-05, 1.6829749E-05, 1.0253250E-05, &
                                     6.4837137E-06, 4.3150590E-06, 3.0681599E-06, 2.3424317E-06, 1.8950741E-06, &
                                     1.5906028E-06, 1.3826408E-06, 1.2976860E-06, 1.5505448E-06, 4.2260924E-06, &
                                     1.8974534E-05, 5.7443782E-05, 1.0313369E-04, 1.2630235E-04, 1.3161768E-04, &
                                     1.4924604E-04, 2.2012185E-04]
      scav_coeff_num(1:ncol, 5, 4) = [ &
                                     9.9093835E-05, 5.4825946E-05, 3.1268211E-05, 1.8348980E-05, 1.1097084E-05, &
                                     6.9638621E-06, 4.5910498E-06, 3.2242133E-06, 2.4262125E-06, 1.9387085E-06, &
                                     1.6203443E-06, 1.4226955E-06, 1.4330785E-06, 2.4009823E-06, 7.7996029E-06, &
                                     2.5431806E-05, 5.9330725E-05, 9.7138366E-05, 1.2153435E-04, 1.3544791E-04, &
                                     1.6558163E-04, 2.6057389E-04]
      scav_coeff_num(1:ncol, 5, 5) = [ &
                                     1.1303342E-04, 6.1687094E-05, 3.4777747E-05, 2.0211323E-05, 1.2118403E-05, &
                                     7.5386978E-06, 4.9188563E-06, 3.4095557E-06, 2.5278463E-06, 1.9947970E-06, &
                                     1.6622435E-06, 1.5055617E-06, 1.8033157E-06, 3.9106121E-06, 1.1598594E-05, &
                                     3.0224368E-05, 6.0461811E-05, 9.3347932E-05, 1.1876960E-04, 1.4115248E-04, &
                                     1.8747535E-04, 3.1630875E-04]
      scav_coeff_num(1:ncol, 6, 1) = [ &
                                     7.9724125E-05, 4.5462621E-05, 2.6569414E-05, 1.5918739E-05, 9.8312923E-06, &
                                     6.3427052E-06, 4.3654193E-06, 3.2669310E-06, 2.6522982E-06, 2.2579003E-06, &
                                     1.9220492E-06, 1.6077184E-06, 1.3988216E-06, 1.4074458E-06, 1.7074232E-06, &
                                     3.9986504E-06, 6.4965897E-05, 1.5758097E-04, 1.6972980E-04, 1.5562152E-04, &
                                     1.5686130E-04, 2.0514049E-04]
      scav_coeff_num(1:ncol, 6, 2) = [ &
                                     8.5349942E-05, 4.8376032E-05, 2.8131064E-05, 1.6782857E-05, 1.0323201E-05, &
                                     6.6286409E-06, 4.5316290E-06, 3.3584494E-06, 2.6940494E-06, 2.2691918E-06, &
                                     1.9263204E-06, 1.6276620E-06, 1.4398123E-06, 1.4672280E-06, 2.3059610E-06, &
                                     1.4290427E-05, 7.0962396E-05, 1.4282305E-04, 1.6478891E-04, 1.5849632E-04, &
                                     1.6578976E-04, 2.2525382E-04]
      scav_coeff_num(1:ncol, 6, 3) = [ &
                                     9.3818301E-05, 5.2713127E-05, 3.0432616E-05, 1.8045069E-05, 1.1035906E-05, &
                                     7.0397556E-06, 4.7690696E-06, 3.4892541E-06, 2.7561133E-06, 2.2913297E-06, &
                                     1.9386878E-06, 1.6562619E-06, 1.5007213E-06, 1.7858313E-06, 5.3708708E-06, &
                                     2.4949868E-05, 7.4467829E-05, 1.3105141E-04, 1.5786268E-04, 1.6191659E-04, &
                                     1.7938623E-04, 2.5734624E-04]
      scav_coeff_num(1:ncol, 6, 4) = [ &
                                     1.0511899E-04, 5.8419720E-05, 3.3421782E-05, 1.9665366E-05, 1.1941218E-05, &
                                     7.5569871E-06, 5.0655865E-06, 3.6529533E-06, 2.8372300E-06, 2.3264001E-06, &
                                     1.9611958E-06, 1.6988872E-06, 1.6792738E-06, 2.9351767E-06, 1.0079747E-05, &
                                     3.3077694E-05, 7.6296483E-05, 1.2309654E-04, 1.5172389E-04, 1.6607541E-04, &
                                     1.9792970E-04, 3.0323612E-04]
      scav_coeff_num(1:ncol, 6, 5) = [ &
                                     1.1957226E-04, 6.5601444E-05, 3.7126813E-05, 2.1646109E-05, 1.3034262E-05, &
                                     8.1746664E-06, 5.4169187E-06, 3.8474987E-06, 2.9374625E-06, 2.3755825E-06, &
                                     1.9986249E-06, 1.7985548E-06, 2.1711692E-06, 4.9318990E-06, 1.4986332E-05, &
                                     3.9020775E-05, 7.7320255E-05, 1.1791794E-04, 1.4785099E-04, 1.7224469E-04, &
                                     2.2270834E-04, 3.6641079E-04]
      scav_coeff_num(1:ncol, 7, 1) = [ &
                                     8.5349386E-05, 4.8801170E-05, 2.8572065E-05, 1.7155635E-05, 1.0649567E-05, &
                                     6.9594621E-06, 4.9177858E-06, 3.8337627E-06, 3.2584061E-06, 2.8725172E-06, &
                                     2.4706447E-06, 2.0180730E-06, 1.6614375E-06, 1.5831273E-06, 1.8789634E-06, &
                                     5.0497533E-06, 8.6796719E-05, 1.9905849E-04, 2.1139246E-04, 1.9220948E-04, &
                                     1.9016876E-04, 2.4168060E-04]
      scav_coeff_num(1:ncol, 7, 2) = [ &
                                     9.1274751E-05, 5.1894385E-05, 3.0241546E-05, 1.8085401E-05, 1.1182516E-05, &
                                     7.2713001E-06, 5.0988747E-06, 3.9299606E-06, 3.2943883E-06, 2.8698488E-06, &
                                     2.4613987E-06, 2.0362388E-06, 1.7138264E-06, 1.6608087E-06, 2.7184571E-06, &
                                     1.8910916E-05, 9.2383956E-05, 1.8080566E-04, 2.0544199E-04, 1.9541800E-04, &
                                     2.0033145E-04, 2.6445765E-04]
      scav_coeff_num(1:ncol, 7, 3) = [ &
                                     1.0016785E-04, 5.6487890E-05, 3.2697233E-05, 1.9441290E-05, 1.1953334E-05, &
                                     7.7185407E-06, 5.3566189E-06, 4.0670274E-06, 3.3493907E-06, 2.8757215E-06, &
                                     2.4588177E-06, 2.0629131E-06, 1.7908008E-06, 2.0950949E-06, 6.8409496E-06, &
                                     3.2545081E-05, 9.5815048E-05, 1.6573371E-04, 1.9689849E-04, 1.9915918E-04, &
                                     2.1575299E-04, 3.0076022E-04]
      scav_coeff_num(1:ncol, 7, 4) = [ &
                                     1.1199026E-04, 6.2512161E-05, 3.5878100E-05, 2.1177951E-05, 1.2930251E-05, &
                                     8.2795754E-06, 5.6772501E-06, 4.2382785E-06, 3.4236941E-06, 2.8952142E-06, &
                                     2.4676256E-06, 2.1065494E-06, 2.0258354E-06, 3.6308759E-06, 1.2983474E-05, &
                                     4.2720770E-05, 9.7511966E-05, 1.5534644E-04, 1.8903482E-04, 2.0358889E-04, &
                                     2.3670731E-04, 3.5261136E-04]
      scav_coeff_num(1:ncol, 7, 5) = [ &
                                     1.2704389E-04, 7.0063927E-05, 3.9807744E-05, 2.3295031E-05, 1.4106594E-05, &
                                     8.9474846E-06, 6.0558363E-06, 4.4418341E-06, 3.5184402E-06, 2.9309969E-06, &
                                     2.4951517E-06, 2.2256544E-06, 2.6739612E-06, 6.2467770E-06, 1.9274675E-05, &
                                     5.0069398E-05, 9.8363493E-05, 1.4842087E-04, 1.8373163E-04, 2.1015577E-04, &
                                     2.6461554E-04, 4.2392017E-04]
      scav_coeff_num(1:ncol, 8, 1) = [ &
                                     9.1692652E-05, 5.2565763E-05, 3.0836030E-05, 1.8566434E-05, 1.1605351E-05, &
                                     7.7145577E-06, 5.6375551E-06, 4.6150183E-06, 4.1277044E-06, 3.7772409E-06, &
                                     3.2907362E-06, 2.6324678E-06, 2.0428613E-06, 1.8150975E-06, 2.0803621E-06, &
                                     6.4140491E-06, 1.1423577E-04, 2.5005939E-04, 2.6284412E-04, 2.3731865E-04, &
                                     2.3066778E-04, 2.8479056E-04]
      scav_coeff_num(1:ncol, 8, 2) = [ &
                                     9.7959062E-05, 5.5862862E-05, 3.2627972E-05, 1.9571493E-05, 1.2186370E-05, &
                                     8.0576191E-06, 5.8368973E-06, 4.7162114E-06, 4.1538019E-06, 3.7520708E-06, &
                                     3.2588204E-06, 2.6455182E-06, 2.1108923E-06, 1.9182494E-06, 3.2365571E-06, &
                                     2.4774448E-05, 1.1907092E-04, 2.2765480E-04, 2.5559325E-04, 2.4084811E-04, &
                                     2.4219796E-04, 3.1049170E-04]
      scav_coeff_num(1:ncol, 8, 3) = [ &
                                     1.0733710E-04, 6.0747602E-05, 3.5258836E-05, 2.1034759E-05, 1.3025039E-05, &
                                     8.5480825E-06, 6.1191638E-06, 4.8595729E-06, 4.1962247E-06, 3.7312988E-06, &
                                     3.2308151E-06, 2.6662951E-06, 2.2092338E-06, 2.5043317E-06, 8.7008471E-06, &
                                     4.2064557E-05, 1.2232563E-04, 2.0857181E-04, 2.4502213E-04, 2.4485513E-04, &
                                     2.5961898E-04, 3.5140014E-04]
      scav_coeff_num(1:ncol, 8, 4) = [ &
                                     1.1975766E-04, 6.7133588E-05, 3.8657828E-05, 2.2904725E-05, 1.4085254E-05, &
                                     9.1610648E-06, 6.4683941E-06, 5.0380189E-06, 4.2575838E-06, 3.7246743E-06, &
                                     3.2158018E-06, 2.7076739E-06, 2.5167911E-06, 4.5317114E-06, 1.6633366E-05, &
                                     5.4745580E-05, 1.2382066E-04, 1.9518872E-04, 2.3498780E-04, 2.4945761E-04, &
                                     2.8318570E-04, 4.0974864E-04]
      scav_coeff_num(1:ncol, 8, 5) = [ &
                                     1.3550362E-04, 7.5108131E-05, 4.2843469E-05, 2.5177929E-05, 1.5358143E-05, &
                                     9.8880078E-06, 6.8787322E-06, 5.2499369E-06, 4.3406813E-06, 3.7372008E-06, &
                                     3.2247186E-06, 2.8480012E-06, 3.3622742E-06, 7.9250703E-06, 2.4642273E-05, &
                                     6.3812480E-05, 1.2443360E-04, 1.8608663E-04, 2.2784961E-04, 2.5630830E-04, &
                                     3.1445892E-04, 4.8989476E-04]
      scav_coeff_num(1:ncol, 9, 1) = [ &
                                     9.8792116E-05, 5.6779523E-05, 3.3377076E-05, 2.0165340E-05, 1.2716829E-05, &
                                     8.6363299E-06, 6.5698626E-06, 5.6782247E-06, 5.3515253E-06, 5.0802428E-06, &
                                     4.4897788E-06, 3.5356410E-06, 2.5934370E-06, 2.1253227E-06, 2.3191693E-06, &
                                     8.1326825E-06, 1.4811269E-04, 3.1242197E-04, 3.2615384E-04, 2.9277995E-04, &
                                     2.7983981E-04, 3.3563328E-04]
      scav_coeff_num(1:ncol, 9, 2) = [ &
                                     1.0544338E-04, 6.0306044E-05, 3.5307202E-05, 2.1256386E-05, 1.3354109E-05, &
                                     9.0171855E-06, 6.7918806E-06, 5.7849874E-06, 5.3623067E-06, 5.0207332E-06, &
                                     4.4219177E-06, 3.5380478E-06, 2.6824747E-06, 2.2642372E-06, 3.8794032E-06, &
                                     3.2061088E-05, 1.5188904E-04, 2.8506476E-04, 3.1722377E-04, 2.9659618E-04, &
                                     2.9287626E-04, 3.6452592E-04]
      scav_coeff_num(1:ncol, 9, 3) = [ &
                                     1.1536992E-04, 6.5518961E-05, 3.8135858E-05, 2.2842187E-05, 1.4271925E-05, &
                                     9.5595652E-06, 7.1041255E-06, 5.9347961E-06, 5.3847496E-06, 4.9588412E-06, &
                                     4.3532545E-06, 3.5464027E-06, 2.8086175E-06, 3.0451035E-06, 1.1013981E-05, &
                                     5.3813812E-05, 1.5489835E-04, 2.6113051E-04, 3.0410095E-04, 3.0077404E-04, &
                                     3.1247105E-04, 4.1044123E-04]
      scav_coeff_num(1:ncol, 9, 4) = [ &
                                     1.2846932E-04, 7.2313373E-05, 4.1781353E-05, 2.4864150E-05, 1.5428915E-05, &
                                     1.0234360E-05, 7.4876356E-06, 6.1199531E-06, 5.4248078E-06, 4.9113047E-06, &
                                     4.2996570E-06, 3.5798139E-06, 3.2072117E-06, 5.6867988E-06, 2.1154166E-05, &
                                     6.9553198E-05, 1.5613987E-04, 2.4408862E-04, 2.9133725E-04, 3.0539580E-04, &
                                     3.3884091E-04, 4.7582175E-04]
      scav_coeff_num(1:ncol, 9, 5) = [ &
                                     1.4500471E-04, 8.0766615E-05, 4.6256664E-05, 2.7315238E-05, 1.6813496E-05, &
                                     1.1030852E-05, 7.9352630E-06, 6.3392103E-06, 5.4877231E-06, 4.8866386E-06, &
                                     4.2770719E-06, 3.7415765E-06, 4.2974550E-06, 1.0043087E-05, 3.1275782E-05, &
                                     8.0721629E-05, 1.5645605E-04, 2.3230325E-04, 2.8185934E-04, 3.1235364E-04, &
                                     3.7369175E-04, 5.6549571E-04]
      scav_coeff_num(1:ncol, 10, 1) = [ &
                                      1.0668249E-04, 6.1463378E-05, 3.6209302E-05, 2.1964856E-05, 1.4000026E-05, &
                                      9.7502489E-06, 7.7563596E-06, 7.0873948E-06, 7.0187546E-06, 6.8895537E-06, &
                                      6.1778736E-06, 4.8166411E-06, 3.3673735E-06, 2.5379705E-06, 2.6037438E-06, &
                                      1.0225711E-05, 1.8917582E-04, 3.8818530E-04, 4.0370676E-04, 3.6073380E-04, &
                                      3.3941444E-04, 3.9553861E-04]
      scav_coeff_num(1:ncol, 10, 2) = [ &
                                      1.1376453E-04, 6.5246178E-05, 3.8294328E-05, 2.3153548E-05, 1.4702876E-05, &
                                      1.0176736E-05, 8.0066347E-06, 7.2007661E-06, 7.0077040E-06, 6.7806597E-06, &
                                      6.0564580E-06, 4.8002906E-06, 3.4836746E-06, 2.7256515E-06, 4.6633233E-06, &
                                      4.0916431E-05, 1.9169289E-04, 3.5490103E-04, 3.9260766E-04, 3.6477089E-04, &
                                      3.5409764E-04, 4.2789103E-04]
      scav_coeff_num(1:ncol, 10, 3) = [ &
                                      1.2430616E-04, 7.0826121E-05, 4.1344811E-05, 2.4878396E-05, 1.5712638E-05, &
                                      1.0781354E-05, 8.3557008E-06, 7.3575511E-06, 7.0011478E-06, 6.6591524E-06, &
                                      5.9268970E-06, 4.7867136E-06, 3.6448548E-06, 3.7500306E-06, 1.3835454E-05, &
                                      6.8077718E-05, 1.9445813E-04, 3.2513199E-04, 3.7626717E-04, 3.6896890E-04, &
                                      3.7603260E-04, 4.7920819E-04]
      scav_coeff_num(1:ncol, 10, 4) = [ &
                                      1.3816896E-04, 7.8078088E-05, 4.5266966E-05, 2.7072683E-05, 1.6981607E-05, &
                                      1.1529595E-05, 8.7805359E-06, 7.5490731E-06, 7.0095688E-06, 6.5517014E-06, &
                                      5.8150636E-06, 4.8035718E-06, 4.1547858E-06, 7.1444065E-06, 2.6662471E-05, &
                                      8.7540722E-05, 1.9543548E-04, 3.0366251E-04, 3.6007800E-04, 3.7338194E-04, &
                                      4.0537294E-04, 5.5213555E-04]
      scav_coeff_num(1:ncol, 10, 5) = [ &
                                      1.5559547E-04, 8.7068839E-05, 5.0067714E-05, 2.9725266E-05, 1.8494879E-05, &
                                      1.2407918E-05, 9.2722022E-06, 7.7745203E-06, 7.0416074E-06, 6.4718497E-06, &
                                      5.7436571E-06, 4.9849042E-06, 5.5425265E-06, 1.2675454E-05, 3.9357848E-05, &
                                      1.0127903E-04, 1.9541929E-04, 2.8860595E-04, 3.4763604E-04, 3.8018721E-04, &
                                      4.4397360E-04, 6.5198786E-04]
      scav_coeff_num(1:ncol, 11, 1) = [ &
                                      1.1539352E-04, 6.6634868E-05, 3.9343653E-05, 2.3973191E-05, 1.5463884E-05, &
                                      1.1069993E-05, 9.2203326E-06, 8.8805685E-06, 9.1855733E-06, 9.2775361E-06, &
                                      8.4328677E-06, 6.5416777E-06, 4.4069240E-06, 3.0727900E-06, 2.9413125E-06, &
                                      1.2683254E-05, 2.3801839E-04, 4.7956381E-04, 4.9822205E-04, 4.4366336E-04, &
                                      4.1140041E-04, 4.6602559E-04]
      scav_coeff_num(1:ncol, 11, 2) = [ &
                                      1.2295401E-04, 7.0701863E-05, 4.1601022E-05, 2.5271844E-05, 1.6242368E-05, &
                                      1.1550867E-05, 9.5054119E-06, 9.0022191E-06, 9.1458000E-06, 9.1022030E-06, &
                                      8.2370802E-06, 6.4961469E-06, 4.5570119E-06, 3.3240749E-06, 5.5976117E-06, &
                                      5.1430260E-05, 2.3927559E-04, 4.3917425E-04, 4.8432479E-04, 4.4780997E-04, &
                                      4.2786747E-04, 5.0210005E-04]
      scav_coeff_num(1:ncol, 11, 3) = [ &
                                      1.3417991E-04, 7.6689226E-05, 4.4898431E-05, 2.7153186E-05, 1.7357909E-05, &
                                      1.2229246E-05, 9.8993285E-06, 9.1671632E-06, 9.1005431E-06, 8.9000127E-06, &
                                      8.0226658E-06, 6.4484987E-06, 4.7605542E-06, 4.6437911E-06, 1.7202965E-05, &
                                      8.5096620E-05, 2.4191831E-04, 4.0243329E-04, 4.6392566E-04, 4.5180174E-04, &
                                      4.5229115E-04, 5.5919518E-04]
      scav_coeff_num(1:ncol, 11, 4) = [ &
                                      1.4889389E-04, 8.4449836E-05, 4.9128745E-05, 2.9541275E-05, 1.8755280E-05, &
                                      1.3063894E-05, 1.0373743E-05, 9.3652155E-06, 9.0660440E-06, 8.7109080E-06, &
                                      7.8293165E-06, 6.4376866E-06, 5.4027459E-06, 8.9406897E-06, 3.3253740E-05, &
                                      1.0907787E-04, 2.4269200E-04, 3.7565957E-04, 4.4344975E-04, 4.5567988E-04, &
                                      4.8473041E-04, 6.4014377E-04]
      scav_coeff_num(1:ncol, 11, 5) = [ &
                                      1.6731714E-04, 9.4039258E-05, 5.4292292E-05, 3.2420295E-05, 2.0415574E-05, &
                                      1.4037694E-05, 1.0917321E-05, 9.5960199E-06, 9.0553942E-06, 8.5552073E-06, &
                                      7.6884774E-06, 6.6348011E-06, 7.1444391E-06, 1.5881557E-05, 4.9051677E-05, &
                                      1.2595532E-04, 2.4235011E-04, 3.5666206E-04, 4.2727856E-04, 4.6196571E-04, &
                                      5.2719166E-04, 7.5075406E-04]
      scav_coeff_num(1:ncol, 12, 1) = [ &
                                      1.2494917E-04, 7.2307336E-05, 4.2786618E-05, 2.6191647E-05, 1.7104993E-05, &
                                      1.2587535E-05, 1.0949735E-05, 1.1043588E-05, 1.1839291E-05, 1.2237397E-05, &
                                      1.1256558E-05, 8.7191548E-06, 5.7214585E-06, 3.7361222E-06, 3.3352825E-06, &
                                      1.5459683E-05, 2.9499600E-04, 5.8890539E-04, 6.1276249E-04, 5.4442619E-04, &
                                      4.9811865E-04, 5.4883058E-04]
      scav_coeff_num(1:ncol, 12, 2) = [ &
                                      1.3303718E-04, 7.6687159E-05, 4.5234092E-05, 2.7612668E-05, 1.7969195E-05, &
                                      1.3131666E-05, 1.1276526E-05, 1.1175861E-05, 1.1764734E-05, 1.1979072E-05, &
                                      1.0965209E-05, 8.6329653E-06, 5.9111077E-06, 4.0657348E-06, 6.6797894E-06, &
                                      6.3616477E-05, 2.9531134E-04, 5.3999983E-04, 5.9526487E-04, 5.4850804E-04, &
                                      5.1649700E-04, 5.8887651E-04]
      scav_coeff_num(1:ncol, 12, 3) = [ &
                                      1.4501849E-04, 8.3123375E-05, 4.8804002E-05, 2.9668119E-05, 1.9204426E-05, &
                                      1.3895527E-05, 1.1723826E-05, 1.1351006E-05, 1.1672051E-05, 1.1675653E-05, &
                                      1.0641543E-05, 8.5379665E-06, 6.1633677E-06, 5.7314359E-06, 2.1127042E-05, &
                                      1.0504138E-04, 2.9813606E-04, 4.9499265E-04, 5.6975328E-04, 5.5196652E-04, &
                                      5.4352825E-04, 6.5209409E-04]
      scav_coeff_num(1:ncol, 12, 4) = [ &
                                      1.6067398E-04, 9.1445044E-05, 5.3374608E-05, 3.2271756E-05, 2.0746774E-05, &
                                      1.4829829E-05, 1.2256670E-05, 1.1556672E-05, 1.1584304E-05, 1.1383610E-05, &
                                      1.0342967E-05, 8.4872951E-06, 6.9578161E-06, 1.1085603E-05, 4.0987899E-05, &
                                      1.3448111E-04, 2.9887620E-04, 4.6193416E-04, 5.4393416E-04, 5.5485718E-04, &
                                      5.7913795E-04, 7.4147685E-04]
      scav_coeff_num(1:ncol, 12, 5) = [ &
                                      1.8020256E-04, 1.0169590E-04, 5.8939114E-05, 3.5402528E-05, 2.2572669E-05, &
                                      1.5913117E-05, 1.2860678E-05, 1.1792902E-05, 1.1520021E-05, 1.1131773E-05, &
                                      1.0111698E-05, 8.6956159E-06, 9.1114355E-06, 1.9688920E-05, 6.0483434E-05, &
                                      1.5518415E-04, 2.9828283E-04, 4.3824814E-04, 5.2310552E-04, 5.6012184E-04, &
                                      6.2548539E-04, 8.6332247E-04]
      scav_coeff_num(1:ncol, 13, 1) = [ &
                                      1.3536731E-04, 7.8489564E-05, 4.6539544E-05, 2.8613090E-05, 1.8904302E-05, &
                                      1.4266675E-05, 1.2885652E-05, 1.3491601E-05, 1.4871935E-05, 1.5650334E-05, &
                                      1.4540484E-05, 1.1271453E-05, 7.2700341E-06, 4.5131877E-06, 3.7828833E-06, &
                                      1.8472238E-05, 3.6012934E-04, 7.1861299E-04, 7.5071694E-04, 6.6626916E-04, &
                                      6.0222583E-04, 6.4593161E-04]
      scav_coeff_num(1:ncol, 13, 2) = [ &
                                      1.4403284E-04, 8.3211160E-05, 4.9194688E-05, 3.0168240E-05, 1.9863322E-05, &
                                      1.4881888E-05, 1.3260462E-05, 1.3637377E-05, 1.4758848E-05, 1.5296439E-05, &
                                      1.4136228E-05, 1.1134224E-05, 7.5030554E-06, 4.9329737E-06, 7.8918994E-06, &
                                      7.7394984E-05, 3.6028370E-04, 6.5952598E-04, 7.2860482E-04, 6.7002698E-04, &
                                      6.2262547E-04, 6.9017947E-04]
      scav_coeff_num(1:ncol, 13, 3) = [ &
                                      1.5684208E-04, 9.0138123E-05, 5.3062429E-05, 3.2414579E-05, 2.1230911E-05, &
                                      1.5741695E-05, 1.3769075E-05, 1.3825672E-05, 1.4613087E-05, 1.4876071E-05, &
                                      1.3683335E-05, 1.0979979E-05, 7.8082223E-06, 6.9886410E-06, 2.5582449E-05, &
                                      1.2798637E-04, 3.6385330E-04, 6.0480933E-04, 6.9667277E-04, 6.7249432E-04, &
                                      6.5233931E-04, 7.5982034E-04]
      scav_coeff_num(1:ncol, 13, 4) = [ &
                                      1.7353108E-04, 9.9073860E-05, 5.8005215E-05, 3.5254602E-05, 2.2933504E-05, &
                                      1.6787657E-05, 1.4368765E-05, 1.4041276E-05, 1.4465178E-05, 1.4464727E-05, &
                                      1.3260146E-05, 1.0878654E-05, 8.7722536E-06, 1.3551026E-05, 4.9875366E-05, &
                                      1.6398331E-04, 3.6488636E-04, 5.6439525E-04, 6.6422828E-04, 6.7378579E-04, &
                                      6.9111357E-04, 8.5796880E-04]
      scav_coeff_num(1:ncol, 13, 5) = [ &
                                      1.9427565E-04, 1.1004966E-04, 6.4008637E-05, 3.8661549E-05, 2.4942314E-05, &
                                      1.7993341E-05, 1.5041527E-05, 1.4284490E-05, 1.4339786E-05, 1.4101128E-05, &
                                      1.2921397E-05, 1.1094823E-05, 1.1394633E-05, 2.4078746E-05, 7.3724120E-05, &
                                      1.8933088E-04, 3.6421464E-04, 5.3520643E-04, 6.3763001E-04, 6.7736206E-04, &
                                      7.4126190E-04, 9.9139497E-04]
      scav_coeff_num(1:ncol, 14, 1) = [ &
                                      1.4665935E-04, 8.5185650E-05, 5.0598777E-05, 3.1222408E-05, 2.0827933E-05, &
                                      1.6044155E-05, 1.4923451E-05, 1.6069650E-05, 1.8079433E-05, 1.9282590E-05, &
                                      1.8060954E-05, 1.4029343E-05, 8.9570951E-06, 5.3658472E-06, 4.2744094E-06, &
                                      2.1604734E-05, 4.3298373E-04, 8.7099626E-04, 9.1572173E-04, 8.1279857E-04, &
                                      7.2670573E-04, 7.5954936E-04]
      scav_coeff_num(1:ncol, 14, 2) = [ &
                                      1.5595289E-04, 9.0277877E-05, 5.3478455E-05, 3.2922068E-05, 2.1888874E-05, &
                                      1.6736023E-05, 1.5350914E-05, 1.6232054E-05, 1.7927734E-05, 1.8828007E-05, &
                                      1.7534766E-05, 1.3834301E-05, 9.2342669E-06, 5.8819119E-06, 9.1995094E-06, &
                                      9.2575177E-05, 4.3438741E-04, 7.9980056E-04, 8.8772650E-04, 8.1586101E-04, &
                                      7.4920928E-04, 8.0820503E-04]
      scav_coeff_num(1:ncol, 14, 3) = [ &
                                      1.6966326E-04, 9.7737356E-05, 5.7668398E-05, 3.5374241E-05, 2.3398910E-05, &
                                      1.7699458E-05, 1.5926910E-05, 1.6437210E-05, 1.7728268E-05, 1.8284150E-05, &
                                      1.6942272E-05, 1.3613214E-05, 9.5933486E-06, 8.3590724E-06, 3.0502614E-05, &
                                      1.5387786E-04, 4.3960973E-04, 7.3381076E-04, 8.4777032E-04, 8.1671283E-04, &
                                      7.8161980E-04, 8.8451501E-04]
      scav_coeff_num(1:ncol, 14, 4) = [ &
                                      1.8747867E-04, 1.0734003E-04, 6.3014133E-05, 3.8469410E-05, 2.5274204E-05, &
                                      1.8866220E-05, 1.6600238E-05, 1.6666449E-05, 1.7519031E-05, 1.7746613E-05, &
                                      1.6384477E-05, 1.3454668E-05, 1.0739940E-05, 1.6266619E-05, 5.9864929E-05, &
                                      1.9769470E-04, 4.4147345E-04, 6.8490878E-04, 8.0716586E-04, 8.1559767E-04, &
                                      8.2345186E-04, 9.9166039E-04]
      scav_coeff_num(1:ncol, 14, 5) = [ &
                                      2.0955089E-04, 1.1910414E-04, 6.9493229E-05, 4.2174785E-05, 2.7480431E-05, &
                                      2.0204529E-05, 1.7348882E-05, 1.6920112E-05, 1.7331041E-05, 1.7264669E-05, &
                                      1.5929752E-05, 1.3678807E-05, 1.3883738E-05, 2.8979219E-05, 8.8771825E-05, &
                                      2.2865011E-04, 4.4103307E-04, 6.4935828E-04, 7.7348573E-04, 8.1662052E-04, &
                                      8.7717702E-04, 1.1368523E-03]
      scav_coeff_num(1:ncol, 15, 1) = [ &
                                      1.5882892E-04, 9.2394684E-05, 5.4956433E-05, 3.3998998E-05, 2.2832537E-05, &
                                      1.7839485E-05, 1.6928826E-05, 1.8576160E-05, 2.1192169E-05, 2.2819966E-05, &
                                      2.1511549E-05, 1.6756139E-05, 1.0645875E-05, 6.2379302E-06, 4.7947191E-06, &
                                      2.4716191E-05, 5.1252354E-04, 1.0480100E-03, 1.1114708E-03, 9.8786049E-04, &
                                      8.7479179E-04, 8.9209396E-04]
      scav_coeff_num(1:ncol, 15, 2) = [ &
                                      1.6880098E-04, 9.7885945E-05, 5.8076427E-05, 3.5851663E-05, 2.3999793E-05, &
                                      1.8610448E-05, 1.7411059E-05, 1.8758192E-05, 2.1006204E-05, 2.2269408E-05, &
                                      2.0866015E-05, 1.6502176E-05, 1.0964374E-05, 6.8485932E-06, 1.0554259E-05, &
                                      1.0884124E-04, 5.1739297E-04, 9.6253946E-04, 1.0760272E-03, 9.8971173E-04, &
                                      8.9944117E-04, 9.4533247E-04]
      scav_coeff_num(1:ncol, 15, 3) = [ &
                                      1.8348569E-04, 1.0591904E-04, 6.2611430E-05, 3.8522065E-05, 2.5658837E-05, &
                                      1.9681283E-05, 1.8057645E-05, 1.8984381E-05, 2.0758593E-05, 2.1607689E-05, &
                                      2.0136958E-05, 1.6213322E-05, 1.1374678E-05, 9.7613540E-06, 3.5777963E-05, &
                                      1.8249535E-04, 5.2561355E-04, 8.8365361E-04, 1.0261126E-03, 9.8811679E-04, &
                                      9.3448027E-04, 1.0284904E-03]
      scav_coeff_num(1:ncol, 15, 4) = [ &
                                      2.0252034E-04, 1.1624069E-04, 6.8389078E-05, 4.1888218E-05, 2.7715403E-05, &
                                      2.0973918E-05, 1.8808799E-05, 1.9232203E-05, 2.0494048E-05, 2.0949673E-05, &
                                      1.9447680E-05, 1.5997472E-05, 1.2711058E-05, 1.9125985E-05, 7.0833831E-05, &
                                      2.3555018E-04, 5.2911775E-04, 8.2512184E-04, 9.7554533E-04, 9.8355336E-04, &
                                      9.7913520E-04, 1.1447442E-03]
      scav_coeff_num(1:ncol, 15, 5) = [ &
                                      2.2603182E-04, 1.2885552E-04, 7.5378608E-05, 4.5911166E-05, 3.0129600E-05, &
                                      2.2451216E-05, 1.9638398E-05, 1.9501771E-05, 2.0249545E-05, 2.0354880E-05, &
                                      1.8880410E-05, 1.6234823E-05, 1.6421328E-05, 3.4268984E-05, 1.0553376E-04, &
                                      2.7322519E-04, 5.2939817E-04, 7.8234491E-04, 9.3326632E-04, 9.8092890E-04, &
                                      1.0360442E-03, 1.3016995E-03]
      scav_coeff_num(1:ncol, 16, 1) = [ &
                                      1.7186835E-04, 1.0010923E-04, 5.9600946E-05, 3.6920208E-05, 2.4873424E-05, &
                                      1.9570507E-05, 1.8764058E-05, 2.0802573E-05, 2.3928196E-05, 2.5930042E-05, &
                                      2.4564003E-05, 1.9194901E-05, 1.2185828E-05, 7.0666815E-06, 5.3265629E-06, &
                                      2.7653628E-05, 5.9695781E-04, 1.2508364E-03, 1.3413472E-03, 1.1952641E-03, &
                                      1.0497638E-03, 1.0460104E-03]
      scav_coeff_num(1:ncol, 16, 2) = [ &
                                      1.8256876E-04, 1.0602710E-04, 6.2975722E-05, 3.8932305E-05, 2.6148484E-05, &
                                      2.0419606E-05, 1.9300318E-05, 2.1006777E-05, 2.3716636E-05, 2.5298478E-05, &
                                      2.3814298E-05, 1.8887488E-05, 1.2539366E-05, 7.7619342E-06, 1.1899437E-05, &
                                      1.2574140E-04, 6.0847266E-04, 1.1487602E-03, 1.2965602E-03, 1.1952081E-03, &
                                      1.0765423E-03, 1.1039675E-03]
      scav_coeff_num(1:ncol, 16, 3) = [ &
                                      1.9830004E-04, 1.1467371E-04, 6.7876796E-05, 4.1830528E-05, 2.7959290E-05, &
                                      2.1597345E-05, 2.0017400E-05, 2.1258246E-05, 2.3432576E-05, 2.4537268E-05, &
                                      2.2966267E-05, 1.8537637E-05, 1.2994640E-05, 1.1103748E-05, 4.1257910E-05, &
                                      2.1340643E-04, 6.2156102E-04, 1.0554041E-03, 1.2344059E-03, 1.1900867E-03, &
                                      1.1140333E-03, 1.1940699E-03]
      scav_coeff_num(1:ncol, 16, 4) = [ &
                                      2.1864553E-04, 1.2576485E-04, 7.4113087E-05, 4.5480246E-05, 3.0201469E-05, &
                                      2.3016439E-05, 2.0847500E-05, 2.1530497E-05, 2.3125963E-05, 2.3777738E-05, &
                                      2.2162856E-05, 1.8271798E-05, 1.4521280E-05, 2.2001113E-05, 8.2579368E-05, &
                                      2.7724084E-04, 6.2784800E-04, 9.8617489E-04, 1.1718121E-03, 1.1807641E-03, &
                                      1.1611168E-03, 1.3194014E-03]
      scav_coeff_num(1:ncol, 16, 5) = [ &
                                      2.4370633E-04, 1.3929101E-04, 8.1645293E-05, 4.9836410E-05, 3.2829856E-05, &
                                      2.4634787E-05, 2.1760622E-05, 2.1823121E-05, 2.2838706E-05, 2.3088318E-05, &
                                      2.1499368E-05, 1.8533009E-05, 1.8831595E-05, 3.9788512E-05, 1.2380637E-04, &
                                      3.2288479E-04, 6.2956374E-04, 9.3536230E-04, 1.1192288E-03, 1.1731456E-03, &
                                      1.2206158E-03, 1.4879005E-03]
      scav_coeff_num(1:ncol, 17, 1) = [ &
                                      1.8574906E-04, 1.0831034E-04, 6.4515524E-05, 3.9963546E-05, 2.6912086E-05, &
                                      2.1169005E-05, 2.0315226E-05, 2.2575256E-05, 2.6050481E-05, 2.8330224E-05, &
                                      2.6935913E-05, 2.1121843E-05, 1.3444200E-05, 7.7960982E-06, 5.8543996E-06, &
                                      3.0267476E-05, 6.8361873E-04, 1.4792784E-03, 1.6077858E-03, 1.4382459E-03, &
                                      1.2545275E-03, 1.2234476E-03]
      scav_coeff_num(1:ncol, 17, 2) = [ &
                                      1.9722585E-04, 1.1468097E-04, 6.8158062E-05, 4.2139555E-05, 2.8293859E-05, &
                                      2.2092274E-05, 2.0902144E-05, 2.2803501E-05, 2.5825369E-05, 2.7641216E-05, &
                                      2.6108162E-05, 2.0772531E-05, 1.3823830E-05, 8.5587295E-06, 1.3176970E-05, &
                                      1.4268626E-04, 7.0600414E-04, 1.3582534E-03, 1.5514234E-03, 1.4353734E-03, &
                                      1.2833372E-03, 1.2862062E-03]
      scav_coeff_num(1:ncol, 17, 3) = [ &
                                      2.1407326E-04, 1.2397895E-04, 7.3444102E-05, 4.5272521E-05, 3.0255772E-05, &
                                      2.3372615E-05, 2.1686563E-05, 2.3084202E-05, 2.5521418E-05, 2.6809602E-05, &
                                      2.5171408E-05, 2.0375528E-05, 1.4315017E-05, 1.2301304E-05, 4.6755695E-05, &
                                      2.4592078E-04, 7.2640836E-04, 1.2490684E-03, 1.4744255E-03, 1.4253645E-03, &
                                      1.3229640E-03, 1.3832461E-03]
      scav_coeff_num(1:ncol, 17, 4) = [ &
                                      2.3581777E-04, 1.3588750E-04, 8.0163166E-05, 4.9215243E-05, 3.2684074E-05, &
                                      2.4914713E-05, 2.2593812E-05, 2.3387101E-05, 2.5191770E-05, 2.5978704E-05, &
                                      2.4283597E-05, 2.0073672E-05, 1.6025033E-05, 2.4759716E-05, 9.4811544E-05, &
                                      3.2213236E-04, 7.3699805E-04, 1.1682699E-03, 1.3975289E-03, 1.4096827E-03, &
                                      1.3718897E-03, 1.5174536E-03]
      scav_coeff_num(1:ncol, 17, 5) = [ &
                                      2.6253385E-04, 1.5038253E-04, 8.8267299E-05, 5.3916816E-05, 3.5528883E-05, &
                                      2.6672232E-05, 2.3590563E-05, 2.3710993E-05, 2.4881737E-05, 2.5223633E-05, &
                                      2.3551357E-05, 2.0373495E-05, 2.0953136E-05, 4.5353892E-05, 1.4325229E-04, &
                                      3.7709579E-04, 7.4112679E-04, 1.1087561E-03, 1.3327998E-03, 1.3954653E-03, &
                                      1.4331547E-03, 1.6970246E-03]
      scav_coeff_num(1:ncol, 18, 1) = [ &
                                      2.0039841E-04, 1.1695429E-04, 6.9671293E-05, 4.3104732E-05, 2.8919381E-05, &
                                      2.2590542E-05, 2.1511092E-05, 2.3785801E-05, 2.7409312E-05, 2.9839122E-05, &
                                      2.8442731E-05, 2.2387971E-05, 1.4330900E-05, 8.3873659E-06, 6.3670326E-06, &
                                      3.2427266E-05, 7.6893481E-04, 1.7309666E-03, 1.9112638E-03, 1.7185348E-03, &
                                      1.4908465E-03, 1.4256391E-03]
      scav_coeff_num(1:ncol, 18, 2) = [ &
                                      2.1269528E-04, 1.2380119E-04, 7.3592635E-05, 4.5447347E-05, 3.0404810E-05, &
                                      2.3581850E-05, 2.2143358E-05, 2.4039194E-05, 2.7184516E-05, 2.9121588E-05, &
                                      2.7570317E-05, 2.2012692E-05, 1.4726258E-05, 9.1953700E-06, 1.4333805E-05, &
                                      1.5896277E-04, 8.0738796E-04, 1.5888883E-03, 1.8408059E-03, 1.7117077E-03, &
                                      1.5214820E-03, 1.4932076E-03]
      scav_coeff_num(1:ncol, 18, 3) = [ &
                                      2.3072210E-04, 1.3378451E-04, 7.9279754E-05, 4.8819472E-05, 3.2514521E-05, &
                                      2.4957701E-05, 2.2989616E-05, 2.4352456E-05, 2.6879946E-05, 2.8255381E-05, &
                                      2.6583482E-05, 2.1587431E-05, 1.5242636E-05, 1.3289646E-05, 5.2055464E-05, &
                                      2.7905428E-04, 8.3811542E-04, 1.4629569E-03, 1.7461338E-03, 1.6951643E-03, &
                                      1.5627590E-03, 1.5970430E-03]
      scav_coeff_num(1:ncol, 18, 4) = [ &
                                      2.5394583E-04, 1.4655349E-04, 8.6502242E-05, 5.3061605E-05, 3.5126394E-05, &
                                      2.6616310E-05, 2.3969910E-05, 2.4691997E-05, 2.6549615E-05, 2.7390304E-05, &
                                      2.5649197E-05, 2.1268091E-05, 1.7121554E-05, 2.7279679E-05, 1.0714773E-04, &
                                      3.6917994E-04, 8.5490872E-04, 1.3700704E-03, 1.6525605E-03, 1.6712508E-03, &
                                      1.6127215E-03, 1.7397159E-03]
      scav_coeff_num(1:ncol, 18, 5) = [ &
                                      2.8241354E-04, 1.6206905E-04, 9.5203458E-05, 5.8117341E-05, 3.8186561E-05, &
                                      2.8508065E-05, 2.5048332E-05, 2.5055574E-05, 2.6240405E-05, 2.6605346E-05, &
                                      2.4882915E-05, 2.1622805E-05, 2.2664712E-05, 5.0767648E-05, 1.6337550E-04, &
                                      4.3483856E-04, 8.6270419E-04, 1.3014473E-03, 1.5738129E-03, 1.6486048E-03, &
                                      1.6746791E-03, 1.9295927E-03]
      scav_coeff_num(1:ncol, 19, 1) = [ &
                                      2.1565206E-04, 1.2594481E-04, 7.5011083E-05, 4.6308676E-05, 3.0871541E-05, &
                                      2.3814869E-05, 2.2328266E-05, 2.4401710E-05, 2.7958895E-05, 3.0397937E-05, &
                                      2.9020409E-05, 2.2937800E-05, 1.4809764E-05, 8.8231511E-06, 6.8577445E-06, &
                                      3.4034502E-05, 8.4856290E-04, 2.0004484E-03, 2.2488616E-03, 2.0349045E-03, &
                                      1.7581057E-03, 1.6518904E-03]
      scav_coeff_num(1:ncol, 19, 2) = [ &
                                      2.2880313E-04, 1.3328599E-04, 7.9219048E-05, 4.8818556E-05, 3.2456094E-05, &
                                      2.4866828E-05, 2.2999444E-05, 2.4680571E-05, 2.7748340E-05, 2.9682253E-05, &
                                      2.8139299E-05, 2.2554499E-05, 1.5210396E-05, 9.6526998E-06, 1.5325956E-05, &
                                      1.7376963E-04, 9.0893417E-04, 1.8358031E-03, 2.1616463E-03, 2.0227819E-03, &
                                      1.7902314E-03, 1.7241690E-03]
      scav_coeff_num(1:ncol, 19, 3) = [ &
                                      2.4805874E-04, 1.4398070E-04, 8.5318697E-05, 5.2431327E-05, 3.4708206E-05, &
                                      2.6329380E-05, 2.3900596E-05, 2.5028817E-05, 2.7462663E-05, 2.8819143E-05, &
                                      2.7144022E-05, 2.2122147E-05, 1.5740800E-05, 1.4031497E-05, 5.6921595E-05, &
                                      3.1151798E-04, 9.5340346E-04, 1.6929133E-03, 2.0464519E-03, 1.9978273E-03, &
                                      1.8324833E-03, 1.8344783E-03]
      scav_coeff_num(1:ncol, 19, 4) = [ &
                                      2.7282423E-04, 1.5764308E-04, 9.3059370E-05, 5.6975603E-05, 3.7498531E-05, &
                                      2.8095986E-05, 2.4948240E-05, 2.5410125E-05, 2.7154415E-05, 2.7959077E-05, &
                                      2.6204085E-05, 2.1805588E-05, 1.7766656E-05, 2.9457124E-05, 1.1911238E-04, &
                                      4.1685736E-04, 9.7860633E-04, 1.5879501E-03, 1.9339383E-03, 1.9636227E-03, &
                                      1.8824510E-03, 1.9849474E-03]
      scav_coeff_num(1:ncol, 19, 5) = [ &
                                      3.0311868E-04, 1.7421877E-04, 1.0237578E-04, 6.2389900E-05, 4.0770022E-05, &
                                      3.0114815E-05, 2.6104743E-05, 2.5821041E-05, 2.6870161E-05, 2.7181926E-05, &
                                      2.5440804E-05, 2.2229841E-05, 2.3897132E-05, 5.5824087E-05, 1.8349985E-04, &
                                      4.9448159E-04, 9.9155929E-04, 1.5101947E-03, 1.8394457E-03, 1.9305921E-03, &
                                      1.9437861E-03, 2.1840186E-03]
      scav_coeff_num(1:ncol, 20, 1) = [ &
                                      2.3117793E-04, 1.3508827E-04, 8.0422520E-05, 4.9512619E-05, 3.2738657E-05, &
                                      2.4836794E-05, 2.2782615E-05, 2.4457380E-05, 2.7747804E-05, 3.0061125E-05, &
                                      2.8717522E-05, 2.2803860E-05, 1.4895159E-05, 9.1051123E-06, 7.3215377E-06, &
                                      3.5029800E-05, 9.1771641E-04, 2.2783716E-03, 2.6125524E-03, 2.3813125E-03, &
                                      2.0516697E-03, 1.8982161E-03]
      scav_coeff_num(1:ncol, 20, 2) = [ &
                                      2.4519801E-04, 1.4293091E-04, 8.4918938E-05, 5.2187248E-05, 3.4416258E-05, &
                                      2.5941320E-05, 2.3485789E-05, 2.4761207E-05, 2.7563854E-05, 2.9375459E-05, &
                                      2.7861779E-05, 2.2430083E-05, 1.5291630E-05, 9.9333383E-06, 1.6119397E-05, &
                                      1.8627595E-04, 1.0058884E-03, 2.0906514E-03, 2.5060595E-03, 2.3624462E-03, &
                                      2.0848075E-03, 1.9749486E-03]
      scav_coeff_num(1:ncol, 20, 3) = [ &
                                      2.6570381E-04, 1.5434738E-04, 9.1433893E-05, 5.6037522E-05, 3.6803062E-05, &
                                      2.7480638E-05, 2.4434073E-05, 2.5145633E-05, 2.7314507E-05, 2.8550436E-05, &
                                      2.6897474E-05, 2.2011237E-05, 1.5825474E-05, 1.4515220E-05, 6.1110834E-05, &
                                      3.4174951E-04, 1.0676047E-03, 1.9315513E-03, 2.3678290E-03, 2.3271255E-03, &
                                      2.1271744E-03, 2.0911713E-03]
      scav_coeff_num(1:ncol, 20, 4) = [ &
                                      2.9203809E-04, 1.6891629E-04, 9.9696520E-05, 6.0880755E-05, 3.9763679E-05, &
                                      2.9345145E-05, 2.5542200E-05, 2.5572283E-05, 2.7048785E-05, 2.7731786E-05, &
                                      2.5990477E-05, 2.1716619E-05, 1.7969257E-05, 3.1207893E-05, 1.3014713E-04, &
                                      4.6312908E-04, 1.1035364E-03, 1.8152148E-03, 2.2345387E-03, 2.2805703E-03, &
                                      2.1759195E-03, 2.2484487E-03]
      scav_coeff_num(1:ncol, 20, 5) = [ &
                                      3.2419160E-04, 1.8656788E-04, 1.0963305E-04, 6.6650855E-05, 4.3238688E-05, &
                                      3.1481821E-05, 2.6771672E-05, 2.6036387E-05, 2.6811205E-05, 2.6997425E-05, &
                                      2.5264825E-05, 2.2220460E-05, 2.4630624E-05, 6.0310713E-05, 2.0276195E-04, &
                                      5.5369184E-04, 1.1232561E-03, 1.7288082E-03, 2.1229777E-03, 2.2352627E-03, &
                                      2.2351268E-03, 2.4552003E-03]
      scav_coeff_num(1:ncol, 21, 1) = [ &
                                      2.4638974E-04, 1.4404260E-04, 8.5706955E-05, 5.2605837E-05, 3.4469203E-05, &
                                      2.5651284E-05, 2.2911719E-05, 2.4031742E-05, 2.6891294E-05, 2.8965362E-05, &
                                      2.7665566E-05, 2.2083678E-05, 1.4637903E-05, 9.2465126E-06, 7.7504340E-06, &
                                      3.5392958E-05, 9.7167576E-04, 2.5511446E-03, 2.9877629E-03, 2.7451491E-03, &
                                      2.3612760E-03, 2.1559896E-03]
      scav_coeff_num(1:ncol, 21, 2) = [ &
                                      2.6125961E-04, 1.5237482E-04, 9.0483154E-05, 5.5437260E-05, 3.6231367E-05, &
                                      2.6799548E-05, 2.3639719E-05, 2.4359071E-05, 2.6743502E-05, 2.8332887E-05, &
                                      2.6863855E-05, 2.1734578E-05, 1.5022481E-05, 1.0053560E-05, 1.6688278E-05, &
                                      1.9569951E-04, 1.0926712E-03, 2.3412415E-03, 2.8600408E-03, 2.7181656E-03, &
                                      2.3948154E-03, 2.2367055E-03]
      scav_coeff_num(1:ncol, 21, 3) = [ &
                                      2.8298758E-04, 1.6449579E-04, 9.7400989E-05, 5.9513955E-05, 3.8741495E-05, &
                                      2.8404341E-05, 2.4626680E-05, 2.4779250E-05, 2.6544048E-05, 2.7574741E-05, &
                                      2.5963648E-05, 2.1346922E-05, 1.5550408E-05, 1.4747501E-05, 6.4388139E-05, &
                                      3.6800367E-04, 1.1747206E-03, 2.1678144E-03, 2.6970793E-03, 2.6707087E-03, &
                                      2.4362956E-03, 2.3579748E-03]
      scav_coeff_num(1:ncol, 21, 4) = [ &
                                      3.1085491E-04, 1.7994935E-04, 1.0617008E-04, 6.4642810E-05, 4.1859375E-05, &
                                      3.0354562E-05, 2.5787122E-05, 2.5252537E-05, 2.6336809E-05, 2.6827427E-05, &
                                      2.5121929E-05, 2.1090005E-05, 1.7776913E-05, 3.2465140E-05, 1.3963823E-04, &
                                      5.0550610E-04, 1.2234962E-03, 2.0416021E-03, 2.5420090E-03, 2.6100407E-03, &
                                      2.4824760E-03, 2.5206912E-03]
      scav_coeff_num(1:ncol, 21, 5) = [ &
                                      3.4482428E-04, 1.9865062E-04, 1.1670852E-04, 7.0753653E-05, 4.5523970E-05, &
                                      3.2596925E-05, 2.7082457E-05, 2.5772823E-05, 2.6162841E-05, 2.6164510E-05, &
                                      2.4462410E-05, 2.1676707E-05, 2.4882164E-05, 6.4010661E-05, 2.2013507E-04, &
                                      6.0944219E-04, 1.2515008E-03, 1.9476079E-03, 2.4127763E-03, 2.5509116E-03, &
                                      2.5379672E-03, 2.7331538E-03]
      scav_coeff_num(1:ncol, 22, 1) = [ &
                                      2.6039575E-04, 1.5228683E-04, 9.0560439E-05, 5.5416246E-05, 3.5978089E-05, &
                                      2.6239495E-05, 2.2755823E-05, 2.3221889E-05, 2.5537079E-05, 2.7290033E-05, &
                                      2.6040320E-05, 2.0909476E-05, 1.4106944E-05, 9.2633958E-06, 8.1290462E-06, &
                                      3.5137250E-05, 1.0064044E-03, 2.8015314E-03, 3.3530810E-03, 3.1065367E-03, &
                                      2.6702933E-03, 2.4113047E-03]
      scav_coeff_num(1:ncol, 22, 2) = [ &
                                      2.7604435E-04, 1.6106764E-04, 9.5591630E-05, 5.8388012E-05, 3.7812396E-05, &
                                      2.7421300E-05, 2.3501040E-05, 2.3570065E-05, 2.5431419E-05, 2.6727295E-05, &
                                      2.5313815E-05, 2.0596502E-05, 1.4473826E-05, 1.0033517E-05, 1.7012244E-05, &
                                      2.0139560E-04, 1.1633780E-03, 2.5719962E-03, 3.2032573E-03, 3.0703877E-03, &
                                      2.7035281E-03, 2.4952641E-03]
      scav_coeff_num(1:ncol, 22, 3) = [ &
                                      2.9889130E-04, 1.7383375E-04, 1.0287660E-04, 6.2667682E-05, 4.0428580E-05, &
                                      2.9077944E-05, 2.4517128E-05, 2.4023487E-05, 2.5290427E-05, 2.6056538E-05, &
                                      2.4502081E-05, 2.0253261E-05, 1.4987955E-05, 1.4744084E-05, 6.6546797E-05, &
                                      3.8851178E-04, 1.2678141E-03, 2.3872883E-03, 3.0152277E-03, 3.0095561E-03, &
                                      2.7430632E-03, 2.6203491E-03]
      scav_coeff_num(1:ncol, 22, 4) = [ &
                                      3.2816048E-04, 1.9009674E-04, 1.1210714E-04, 6.8053065E-05, 4.3682932E-05, &
                                      3.1098244E-05, 2.5719750E-05, 2.4541712E-05, 2.5151756E-05, 2.5401769E-05, &
                                      2.3749277E-05, 2.0044944E-05, 1.7257137E-05, 3.3177060E-05, 1.4696676E-04, &
                                      5.4122571E-04, 1.3309402E-03, 2.2534896E-03, 2.8386298E-03, 2.9336818E-03, &
                                      2.7853536E-03, 2.7867097E-03]
      scav_coeff_num(1:ncol, 22, 5) = [ &
                                      3.6378751E-04, 2.0975741E-04, 1.2319388E-04, 7.4470560E-05, 4.7513526E-05, &
                                      3.3429452E-05, 2.7070889E-05, 2.5117344E-05, 2.5052013E-05, 2.4830415E-05, &
                                      2.3176014E-05, 2.0709996E-05, 2.4689649E-05, 6.6712154E-05, 2.3450105E-04, &
                                      6.5818696E-04, 1.3683810E-03, 2.1535597E-03, 2.6921600E-03, 2.8598736E-03, &
                                      2.8356157E-03, 3.0024301E-03]

      scav_coeff_mass(1:ncol, 1, 1) = [ &
                                      5.3745231E-05, 3.0471292E-05, 1.7799763E-05, 1.0682883E-05, 6.5863179E-06, &
                                      4.1856180E-06, 2.7625087E-06, 1.9144385E-06, 1.4094690E-06, 1.1117903E-06, &
                                      9.4386013E-07, 8.6775625E-07, 8.7377425E-07, 9.7328374E-07, 1.2046318E-06, &
                                      2.0995472E-06, 1.8703475E-05, 4.8626027E-05, 5.5368202E-05, 5.4622025E-05, &
                                      6.4438542E-05, 1.0226403E-04]
      scav_coeff_mass(1:ncol, 1, 2) = [ &
                                      4.2612415E-05, 2.4375993E-05, 1.4360801E-05, 8.6965261E-06, 5.4194431E-06, &
                                      3.4929053E-06, 2.3494812E-06, 1.6686298E-06, 1.2649820E-06, 1.0313155E-06, &
                                      9.0977269E-07, 8.7656299E-07, 9.3306178E-07, 1.1239419E-06, 2.3087080E-06, &
                                      1.1354039E-05, 3.4115185E-05, 5.1223839E-05, 5.5116574E-05, 6.1034258E-05, &
                                      8.7268862E-05, 1.6307426E-04]
      scav_coeff_mass(1:ncol, 1, 3) = [ &
                                      3.1305905E-05, 1.8125496E-05, 1.0807173E-05, 6.6326575E-06, 4.2030985E-06, &
                                      2.7701101E-06, 1.9191745E-06, 1.4142874E-06, 1.1196859E-06, 9.6005398E-07, &
                                      8.9997776E-07, 9.3756281E-07, 1.2074771E-06, 2.9182819E-06, 1.0455557E-05, &
                                      2.7065133E-05, 4.4472562E-05, 5.3813904E-05, 6.1775700E-05, 8.6073684E-05, &
                                      1.5637231E-04, 3.3835058E-04]
      scav_coeff_mass(1:ncol, 1, 4) = [ &
                                      2.2049785E-05, 1.2951012E-05, 7.8401186E-06, 4.8995734E-06, 3.1787770E-06, &
                                      2.1616183E-06, 1.5591140E-06, 1.2068049E-06, 1.0129911E-06, 9.3394044E-07, &
                                      9.8784714E-07, 1.4688329E-06, 3.7425976E-06, 1.0822082E-05, 2.4302812E-05, &
                                      3.9804694E-05, 5.1701430E-05, 6.3501816E-05, 9.1095100E-05, 1.6668155E-04, &
                                      3.6185721E-04, 8.5068497E-04]
      scav_coeff_mass(1:ncol, 1, 5) = [ &
                                      1.5205923E-05, 9.0816803E-06, 5.6033221E-06, 3.5866887E-06, 2.4019250E-06, &
                                      1.7025239E-06, 1.2939047E-06, 1.0687080E-06, 9.8299404E-07, 1.1063105E-06, &
                                      1.8864481E-06, 4.6635416E-06, 1.1540412E-05, 2.3269516E-05, 3.7223857E-05, &
                                      5.0376600E-05, 6.6239498E-05, 1.0053311E-04, 1.8959552E-04, 4.1671403E-04, &
                                      9.7631825E-04, 2.2715296E-03]
      scav_coeff_mass(1:ncol, 2, 1) = [ &
                                      5.6063579E-05, 3.1847880E-05, 1.8617014E-05, 1.1173214E-05, 6.8888188E-06, &
                                      4.3829845E-06, 2.9033928E-06, 2.0262403E-06, 1.5054641E-06, 1.1957061E-06, &
                                      1.0147457E-06, 9.2645788E-07, 9.2785140E-07, 1.0367287E-06, 1.2969745E-06, &
                                      2.5796045E-06, 2.6417692E-05, 6.3026124E-05, 6.9206335E-05, 6.7067651E-05, &
                                      7.7368874E-05, 1.2006288E-04]
      scav_coeff_mass(1:ncol, 2, 2) = [ &
                                      4.4471208E-05, 2.5480453E-05, 1.5019381E-05, 9.0959547E-06, 5.6712865E-06, &
                                      3.6632230E-06, 2.4764776E-06, 1.7727172E-06, 1.3547692E-06, 1.1084590E-06, &
                                      9.7460454E-07, 9.3355668E-07, 9.9316647E-07, 1.2163740E-06, 2.9000762E-06, &
                                      1.5607006E-05, 4.5154477E-05, 6.5078995E-05, 6.8267674E-05, 7.3942370E-05, &
                                      1.0321728E-04, 1.8971105E-04]
      scav_coeff_mass(1:ncol, 2, 3) = [ &
                                      3.2687488E-05, 1.8948913E-05, 1.1303007E-05, 6.9398664E-06, 4.4043401E-06, &
                                      2.9136490E-06, 2.0317622E-06, 1.5088221E-06, 1.2007138E-06, 1.0291133E-06, &
                                      9.6079036E-07, 1.0019051E-06, 1.3532496E-06, 3.7446933E-06, 1.4059691E-05, &
                                      3.5717509E-05, 5.6969283E-05, 6.7119231E-05, 7.5164192E-05, 1.0214278E-04, &
                                      1.8225454E-04, 3.9112904E-04]
      scav_coeff_mass(1:ncol, 2, 4) = [ &
                                      2.3032570E-05, 1.3541404E-05, 8.2024565E-06, 5.1322591E-06, 3.3396168E-06, &
                                      2.2832421E-06, 1.6582995E-06, 1.2908383E-06, 1.0851219E-06, 9.9924430E-07, &
                                      1.0711032E-06, 1.7253685E-06, 4.8472812E-06, 1.4358597E-05, 3.1841397E-05, &
                                      5.1002928E-05, 6.4651509E-05, 7.7367913E-05, 1.0814713E-04, 1.9431944E-04, &
                                      4.1837226E-04, 9.8078995E-04]
      scav_coeff_mass(1:ncol, 2, 5) = [ &
                                      1.5889885E-05, 9.4993086E-06, 5.8680113E-06, 3.7654846E-06, 2.5331363E-06, &
                                      1.8066192E-06, 1.3808048E-06, 1.1438978E-06, 1.0571559E-06, 1.2323535E-06, &
                                      2.2973392E-06, 6.0542218E-06, 1.5169707E-05, 3.0273798E-05, 4.7544866E-05, &
                                      6.2907743E-05, 8.0573298E-05, 1.1916107E-04, 2.2085457E-04, 4.8168287E-04, &
                                      1.1255995E-03, 2.6169100E-03]
      scav_coeff_mass(1:ncol, 3, 1) = [ &
                                      5.8865063E-05, 3.3511807E-05, 1.9608103E-05, 1.1772057E-05, 7.2633920E-06, &
                                      4.6337085E-06, 3.0897236E-06, 2.1814105E-06, 1.6443404E-06, 1.3198936E-06, &
                                      1.1191012E-06, 1.0090169E-06, 9.9761086E-07, 1.1123232E-06, 1.4030193E-06, &
                                      3.2742033E-06, 3.6751432E-05, 8.1008027E-05, 8.6383969E-05, 8.2415738E-05, &
                                      9.2956540E-05, 1.4082319E-04]
      scav_coeff_mass(1:ncol, 3, 2) = [ &
                                      4.6718202E-05, 2.6817213E-05, 1.5820184E-05, 9.5863739E-06, 5.9863461E-06, &
                                      3.8832942E-06, 2.6480621E-06, 1.9198921E-06, 1.4858752E-06, 1.2221200E-06, &
                                      1.0678705E-06, 1.0103460E-06, 1.0674572E-06, 1.3287404E-06, 3.7143753E-06, &
                                      2.1296290E-05, 5.9238788E-05, 8.2347562E-05, 8.4531440E-05, 8.9650033E-05, &
                                      1.2203916E-04, 2.2020836E-04]
      scav_coeff_mass(1:ncol, 3, 3) = [ &
                                      3.4359274E-05, 1.9948103E-05, 1.1909119E-05, 7.3210088E-06, 4.6607807E-06, &
                                      3.1040085E-06, 2.1879692E-06, 1.6448437E-06, 1.3192211E-06, 1.1288584E-06, &
                                      1.0443371E-06, 1.0841229E-06, 1.5427714E-06, 4.8564869E-06, 1.8809883E-05, &
                                      4.6806845E-05, 7.2653944E-05, 8.3604141E-05, 9.1491523E-05, 1.2118555E-04, &
                                      2.1198388E-04, 4.5059897E-04]
      scav_coeff_mass(1:ncol, 3, 4) = [ &
                                      2.4224351E-05, 1.4261328E-05, 8.6497348E-06, 5.4261269E-06, 3.5501215E-06, &
                                      2.4494350E-06, 1.7990594E-06, 1.4125482E-06, 1.1889538E-06, 1.0895776E-06, &
                                      1.1806136E-06, 2.0650154E-06, 6.3110116E-06, 1.8969859E-05, 4.1486846E-05, &
                                      6.5095674E-05, 8.0719393E-05, 9.4265876E-05, 1.2833979E-04, 2.2605569E-04, &
                                      4.8204648E-04, 1.1262084E-03]
      scav_coeff_mass(1:ncol, 3, 5) = [ &
                                      1.6722850E-05, 1.0013141E-05, 6.2002053E-06, 3.9971731E-06, 2.7101916E-06, &
                                      1.9524972E-06, 1.5053662E-06, 1.2514240E-06, 1.1597288E-06, 1.4001977E-06, &
                                      2.8403835E-06, 7.8777049E-06, 1.9866494E-05, 3.9209357E-05, 6.0530731E-05, &
                                      7.8447370E-05, 9.7994258E-05, 1.4113923E-04, 2.5665878E-04, 5.5480118E-04, &
                                      1.2923894E-03, 3.0018139E-03]
      scav_coeff_mass(1:ncol, 4, 1) = [ &
                                      6.2174955E-05, 3.5477932E-05, 2.0782433E-05, 1.2486514E-05, 7.7171032E-06, &
                                      4.9466042E-06, 3.3333879E-06, 2.3954760E-06, 1.8445944E-06, 1.5035616E-06, &
                                      1.2734719E-06, 1.1266545E-06, 1.0888987E-06, 1.2026400E-06, 1.5243501E-06, &
                                      4.2706047E-06, 5.0354259E-05, 1.0339176E-04, 1.0771350E-04, 1.0135387E-04, &
                                      1.1177062E-04, 1.6504866E-04]
      scav_coeff_mass(1:ncol, 4, 2) = [ &
                                      4.9373664E-05, 2.8398545E-05, 1.6771585E-05, 1.0174995E-05, 6.3726958E-06, &
                                      4.1635909E-06, 2.8780101E-06, 2.1271299E-06, 1.6769869E-06, 1.3899341E-06, &
                                      1.2032660E-06, 1.1154200E-06, 1.1602290E-06, 1.4659631E-06, 4.8225208E-06, &
                                      2.8804503E-05, 7.7081977E-05, 1.0383739E-04, 1.0465185E-04, 1.0878505E-04, &
                                      1.4427413E-04, 2.5509006E-04]
      scav_coeff_mass(1:ncol, 4, 3) = [ &
                                      3.6336622E-05, 2.1133059E-05, 1.2633391E-05, 7.7841450E-06, 4.9823298E-06, &
                                      3.3539922E-06, 2.4036357E-06, 1.8401584E-06, 1.4927363E-06, 1.2739262E-06, &
                                      1.1607486E-06, 1.1905484E-06, 1.7889954E-06, 6.3352703E-06, 2.5004438E-05, &
                                      6.0922175E-05, 9.2280461E-05, 1.0402064E-04, 1.1142176E-04, 1.4377745E-04, &
                                      2.4610724E-04, 5.1740537E-04]
      scav_coeff_mass(1:ncol, 4, 4) = [ &
                                      2.5636893E-05, 1.5119589E-05, 9.1903592E-06, 5.7909781E-06, 3.8226000E-06, &
                                      2.6752631E-06, 1.9983790E-06, 1.5889699E-06, 1.3393063E-06, 1.2160884E-06, &
                                      1.3258378E-06, 2.5119789E-06, 8.2310650E-06, 2.4930036E-05, 5.3752433E-05, &
                                      8.2771822E-05, 1.0063685E-04, 1.1487383E-04, 1.5227678E-04, 2.6246917E-04, &
                                      5.5356416E-04, 1.2880232E-03]
      scav_coeff_mass(1:ncol, 4, 5) = [ &
                                      1.7714692E-05, 1.0632124E-05, 6.6098416E-06, 4.2938386E-06, 2.9476229E-06, &
                                      2.1564451E-06, 1.6840448E-06, 1.4060173E-06, 1.3029455E-06, 1.6240283E-06, &
                                      3.5527390E-06, 1.0247752E-05, 2.5899479E-05, 5.0546911E-05, 7.6816918E-05, &
                                      9.7696509E-05, 1.1918002E-04, 1.6709231E-04, 2.9762442E-04, 6.3681881E-04, &
                                      1.4779011E-03, 3.4286285E-03]
      scav_coeff_mass(1:ncol, 5, 1) = [ &
                                      6.6022410E-05, 3.7763585E-05, 2.2151225E-05, 1.3325469E-05, 8.2593716E-06, &
                                      5.3340584E-06, 3.6516762E-06, 2.6915668E-06, 2.1343235E-06, 1.7764272E-06, &
                                      1.5039112E-06, 1.2973164E-06, 1.2111032E-06, 1.3116243E-06, 1.6630317E-06, &
                                      5.6779651E-06, 6.7986717E-05, 1.3118650E-04, 1.3419538E-04, 1.2472807E-04, &
                                      1.3450674E-04, 1.9334211E-04]
      scav_coeff_mass(1:ncol, 5, 2) = [ &
                                      5.2461067E-05, 3.0238894E-05, 1.7883799E-05, 1.0871149E-05, 6.8414381E-06, &
                                      4.5191418E-06, 3.1866488E-06, 2.4200056E-06, 1.9567670E-06, 1.6394102E-06, &
                                      1.4024814E-06, 1.2623565E-06, 1.2782647E-06, 1.6346185E-06, 6.3102336E-06, &
                                      3.8586728E-05, 9.9546866E-05, 1.3054136E-04, 1.2954456E-04, 1.3211714E-04, &
                                      1.7057435E-04, 2.9496663E-04]
      scav_coeff_mass(1:ncol, 5, 3) = [ &
                                      3.8637482E-05, 2.2515805E-05, 1.3485837E-05, 8.3401847E-06, 5.3830705E-06, &
                                      3.6823789E-06, 2.7025122E-06, 2.1219429E-06, 1.7484391E-06, 1.4873349E-06, &
                                      1.3261173E-06, 1.3310312E-06, 2.1085397E-06, 8.2791366E-06, 3.3000309E-05, &
                                      7.8774565E-05, 1.1677011E-04, 1.2929290E-04, 1.3576924E-04, 1.7061558E-04, &
                                      2.8526600E-04, 5.9226060E-04]
      scav_coeff_mass(1:ncol, 5, 4) = [ &
                                      2.7284249E-05, 1.6127239E-05, 9.8355802E-06, 6.2406433E-06, 4.1750323E-06, &
                                      2.9832494E-06, 2.2821031E-06, 1.8464167E-06, 1.5593669E-06, 1.3963239E-06, &
                                      1.5207807E-06, 3.0966640E-06, 1.0724527E-05, 3.2568514E-05, 6.9258715E-05, &
                                      1.0487119E-04, 1.2530086E-04, 1.4002222E-04, 1.8068889E-04, 3.0423832E-04, &
                                      6.3367996E-04, 1.4673464E-03]
      scav_coeff_mass(1:ncol, 5, 5) = [ &
                                      1.8877709E-05, 1.1368131E-05, 7.1108642E-06, 4.6730790E-06, 3.2671104E-06, &
                                      2.4431948E-06, 1.9421881E-06, 1.6306157E-06, 1.5057215E-06, 1.9237466E-06, &
                                      4.4805973E-06, 1.3301649E-05, 3.3592505E-05, 6.4856458E-05, 9.7177514E-05, &
                                      1.2151264E-04, 1.4495849E-04, 1.9777419E-04, 3.4446996E-04, 7.2855628E-04, &
                                      1.6833721E-03, 3.8996803E-03]
      scav_coeff_mass(1:ncol, 6, 1) = [ &
                                      7.0439935E-05, 4.0388284E-05, 2.3727442E-05, 1.4299750E-05, 8.9024966E-06, &
                                      5.8131340E-06, 4.0692331E-06, 3.1034519E-06, 2.5553863E-06, 2.1836199E-06, &
                                      1.8507485E-06, 1.5492847E-06, 1.3791757E-06, 1.4453949E-06, 1.8218224E-06, &
                                      7.6242244E-06, 9.0523819E-05, 1.6561720E-04, 1.6704913E-04, 1.5357086E-04, &
                                      1.6200857E-04, 2.2641857E-04]
      scav_coeff_mass(1:ncol, 6, 2) = [ &
                                      5.6006685E-05, 3.2354708E-05, 1.9168931E-05, 1.1686650E-05, 7.4075464E-06, &
                                      4.9710926E-06, 3.6033820E-06, 2.8363024E-06, 2.3683610E-06, 2.0127028E-06, &
                                      1.6993240E-06, 1.4726046E-06, 1.4322516E-06, 1.8433755E-06, 8.2766856E-06, &
                                      5.1170271E-05, 1.2766158E-04, 1.6366547E-04, 1.6032605E-04, 1.6058233E-04, &
                                      2.0172137E-04, 3.4054102E-04]
      scav_coeff_mass(1:ncol, 6, 3) = [ &
                                      4.1282153E-05, 2.4110392E-05, 1.4478891E-05, 9.0036195E-06, 5.8826422E-06, &
                                      4.1161813E-06, 3.1195249E-06, 2.5309030E-06, 2.1277281E-06, 1.8046937E-06, &
                                      1.5656560E-06, 1.5208025E-06, 2.5223494E-06, 1.0802177E-05, 4.3214389E-05, &
                                      1.0120838E-04, 1.4723480E-04, 1.6054509E-04, 1.6552390E-04, 2.0253736E-04, &
                                      3.3020389E-04, 6.7593268E-04]
      scav_coeff_mass(1:ncol, 6, 4) = [ &
                                      2.9182718E-05, 1.7297846E-05, 1.0600186E-05, 6.7942896E-06, 4.6331920E-06, &
                                      3.4064296E-06, 2.6888317E-06, 2.2248140E-06, 1.8848193E-06, 1.6575214E-06, &
                                      1.7861090E-06, 3.8566024E-06, 1.3928607E-05, 4.2272406E-05, 8.8743474E-05, &
                                      1.3240423E-04, 1.5579902E-04, 1.7071865E-04, 2.1445291E-04, 3.5214916E-04, &
                                      7.2320599E-04, 1.6652607E-03]
      scav_coeff_mass(1:ncol, 6, 5) = [ &
                                      2.0226878E-05, 1.2236652E-05, 7.7225089E-06, 5.1600602E-06, 3.7004028E-06, &
                                      2.8496167E-06, 2.3181542E-06, 1.9603547E-06, 1.7969520E-06, 2.3272234E-06, &
                                      5.6804210E-06, 1.7201303E-05, 4.3328050E-05, 8.2816938E-05, 1.2254287E-04, &
                                      1.5093285E-04, 1.7633129E-04, 2.3408617E-04, 3.9802276E-04, 8.3088771E-04, &
                                      1.9099949E-03, 4.4170477E-03]
      scav_coeff_mass(1:ncol, 7, 1) = [ &
                                      7.5462290E-05, 4.3373111E-05, 2.5525473E-05, 1.5422099E-05, 9.6620052E-06, &
                                      6.4065020E-06, 4.6198849E-06, 3.6785871E-06, 3.1678175E-06, 2.7911550E-06, &
                                      2.3741637E-06, 1.9255898E-06, 1.6161996E-06, 1.6132992E-06, 2.0044454E-06, &
                                      1.0247291E-05, 1.1894579E-04, 2.0814418E-04, 2.0774270E-04, 1.8912905E-04, &
                                      1.9528823E-04, 2.6511492E-04]
      scav_coeff_mass(1:ncol, 7, 2) = [ &
                                      6.0038726E-05, 3.4763963E-05, 2.0640814E-05, 1.2635970E-05, 8.0905340E-06, &
                                      5.5481249E-06, 4.1691664E-06, 3.4297652E-06, 2.9743206E-06, 2.5720687E-06, &
                                      2.1446289E-06, 1.7789285E-06, 1.6385921E-06, 2.1034137E-06, 1.0830137E-05, &
                                      6.7144536E-05, 1.6262634E-04, 2.0465266E-04, 1.9834164E-04, 1.9530671E-04, &
                                      2.3864055E-04, 3.9260744E-04]
      scav_coeff_mass(1:ncol, 7, 3) = [ &
                                      4.4292677E-05, 2.5932650E-05, 1.5627507E-05, 9.7930809E-06, 6.5074794E-06, &
                                      4.6928411E-06, 3.7041413E-06, 3.1257684E-06, 2.6913716E-06, 2.2791354E-06, &
                                      1.9175302E-06, 1.7827837E-06, 3.0561577E-06, 1.4031368E-05, 5.6117708E-05, &
                                      1.2920370E-04, 1.8499432E-04, 1.9912670E-04, 2.0187539E-04, 2.4053675E-04, &
                                      3.8176829E-04, 7.6921620E-04]
      scav_coeff_mass(1:ncol, 7, 4) = [ &
                                      3.1350519E-05, 1.8647538E-05, 1.1503017E-05, 7.4775889E-06, 5.2327082E-06, &
                                      3.9914931E-06, 3.2741103E-06, 2.7825510E-06, 2.3686302E-06, 2.0405424E-06, &
                                      2.1516604E-06, 4.8369580E-06, 1.7998390E-05, 5.4483254E-05, 1.1306403E-04, &
                                      1.6656545E-04, 1.9343153E-04, 2.0817060E-04, 2.5460779E-04, 4.0709608E-04, &
                                      8.2298034E-04, 1.8827134E-03]
      scav_coeff_mass(1:ncol, 7, 5) = [ &
                                      2.1779887E-05, 1.3257300E-05, 8.4704528E-06, 5.7895173E-06, 4.2923221E-06, &
                                      3.4286882E-06, 2.8679014E-06, 2.4471583E-06, 2.2193904E-06, 2.8728252E-06, &
                                      7.2195375E-06, 2.2131767E-05, 5.5545634E-05, 1.0521923E-04, 1.5401117E-04, &
                                      1.8719289E-04, 2.1449446E-04, 2.7709150E-04, 4.5921707E-04, 9.4470199E-04, &
                                      2.1587925E-03, 4.9822431E-03]
      scav_coeff_mass(1:ncol, 8, 1) = [ &
                                      8.1124862E-05, 4.6739721E-05, 2.7560506E-05, 1.6706750E-05, 1.0556382E-05, &
                                      7.1424042E-06, 5.3470405E-06, 4.4792670E-06, 4.0520341E-06, 3.6892288E-06, &
                                      3.1580461E-06, 2.4874210E-06, 1.9555619E-06, 1.8288691E-06, 2.2158268E-06, &
                                      1.3679779E-05, 1.5431208E-04, 2.6047316E-04, 2.5801885E-04, 2.3289077E-04, &
                                      2.3554581E-04, 3.1039640E-04]
      scav_coeff_mass(1:ncol, 8, 2) = [ &
                                      6.4586011E-05, 3.7485349E-05, 2.2314482E-05, 1.3735902E-05, 8.9143148E-06, &
                                      6.2866636E-06, 4.9373242E-06, 4.2718107E-06, 3.8593594E-06, 3.4034053E-06, &
                                      2.8098296E-06, 2.2281511E-06, 1.9209552E-06, 2.4285659E-06, 1.4080018E-05, &
                                      8.7138758E-05, 2.0580781E-04, 2.5520029E-04, 2.4519144E-04, 2.3762899E-04, &
                                      2.8241322E-04, 4.5204416E-04]
      scav_coeff_mass(1:ncol, 8, 3) = [ &
                                      4.7691832E-05, 2.7999546E-05, 1.6948751E-05, 1.0731157E-05, 7.2909505E-06, &
                                      5.4609028E-06, 4.5218383E-06, 3.9857103E-06, 3.5227172E-06, 2.9847551E-06, &
                                      2.4357788E-06, 2.1494364E-06, 3.7402784E-06, 1.8100176E-05, 7.2220838E-05, &
                                      1.6386726E-04, 2.3158551E-04, 2.4663419E-04, 2.4623546E-04, 2.8577825E-04, &
                                      4.4090589E-04, 8.7288901E-04]
      scav_coeff_mass(1:ncol, 8, 4) = [ &
                                      3.3807059E-05, 2.0194561E-05, 1.2566763E-05, 8.3228389E-06, 6.0196915E-06, &
                                      4.8001355E-06, 4.1126405E-06, 3.5994301E-06, 3.0842909E-06, 2.6027591E-06, &
                                      2.6582676E-06, 6.0899368E-06, 2.3101242E-05, 6.9686434E-05, 1.4318976E-04, &
                                      2.0873848E-04, 2.3972812E-04, 2.5380571E-04, 3.0236797E-04, 4.7007658E-04, &
                                      9.3382004E-04, 2.1203718E-03]
      scav_coeff_mass(1:ncol, 8, 5) = [ &
                                      2.3556632E-05, 1.4453598E-05, 9.3869416E-06, 6.6063761E-06, 5.1020622E-06, &
                                      4.2515719E-06, 3.6677387E-06, 3.1627799E-06, 2.8323825E-06, 3.6110112E-06, &
                                      9.1752366E-06, 2.8296294E-05, 7.0733733E-05, 1.3296059E-04, 1.9285204E-04, &
                                      2.3174067E-04, 2.6085556E-04, 3.2802632E-04, 5.2908475E-04, 1.0708457E-03, &
                                      2.4304500E-03, 5.5957880E-03]
      scav_coeff_mass(1:ncol, 9, 1) = [ &
                                      8.7461745E-05, 5.0509088E-05, 2.9847508E-05, 1.8168284E-05, 1.1605462E-05, &
                                      8.0522544E-06, 6.3002945E-06, 5.5782757E-06, 5.3039705E-06, 4.9876748E-06, &
                                      4.3066891E-06, 3.3124368E-06, 2.4404438E-06, 2.1097768E-06, 2.4620617E-06, &
                                      1.8028442E-05, 1.9771836E-04, 3.2455599E-04, 3.1991924E-04, 2.8661284E-04, &
                                      2.8418951E-04, 3.6336325E-04]
      scav_coeff_mass(1:ncol, 9, 2) = [ &
                                      6.9676370E-05, 4.0537116E-05, 2.4205063E-05, 1.5004147E-05, 9.9051666E-06, &
                                      7.2279543E-06, 5.9696656E-06, 5.4469510E-06, 5.1257107E-06, 4.6123910E-06, &
                                      3.7844900E-06, 2.8800428E-06, 2.3099690E-06, 2.8346282E-06, 1.8125757E-05, &
                                      1.1178826E-04, 2.5872087E-04, 3.1727190E-04, 3.0275522E-04, 2.8912364E-04, &
                                      3.3428923E-04, 5.1980553E-04]
      scav_coeff_mass(1:ncol, 9, 3) = [ &
                                      5.1501781E-05, 3.0328019E-05, 1.8460446E-05, 1.1842503E-05, 8.2706841E-06, &
                                      6.4764502E-06, 5.6498169E-06, 5.2058362E-06, 4.7236944E-06, 4.0137413E-06, &
                                      3.1887309E-06, 2.6620306E-06, 4.6078234E-06, 2.3138887E-05, 9.2050973E-05, &
                                      2.0641237E-04, 2.8876364E-04, 3.0492965E-04, 3.0026009E-04, 3.3961106E-04, &
                                      5.0865769E-04, 9.8766574E-04]
      scav_coeff_mass(1:ncol, 9, 4) = [ &
                                      3.6571637E-05, 2.1957870E-05, 1.3816114E-05, 9.3664117E-06, 7.0473688E-06, &
                                      5.9050092E-06, 5.2939667E-06, 4.7727084E-06, 4.1227798E-06, 3.4161418E-06, &
                                      3.3564581E-06, 7.6718989E-06, 2.9407596E-05, 8.8393537E-05, 1.8018522E-04, &
                                      2.6049314E-04, 2.9646081E-04, 3.0929057E-04, 3.5913648E-04, 5.4218507E-04, &
                                      1.0564707E-03, 2.3784686E-03]
      scav_coeff_mass(1:ncol, 9, 5) = [ &
                                      2.5577743E-05, 1.5851111E-05, 1.0508292E-05, 7.6625165E-06, 6.1993265E-06, &
                                      5.4034901E-06, 4.8104695E-06, 4.1957082E-06, 3.7096793E-06, 4.6023956E-06, &
                                      1.1630948E-05, 3.5907415E-05, 8.9415276E-05, 1.6703088E-04, 2.4050260E-04, &
                                      2.8624490E-04, 3.1704878E-04, 3.8830970E-04, 6.0874462E-04, 1.2100621E-03, &
                                      2.7251344E-03, 6.2567626E-03]
      scav_coeff_mass(1:ncol, 10, 1) = [ &
                                       9.4503917E-05, 5.4700162E-05, 3.2399813E-05, 1.9819515E-05, 1.2826775E-05, &
                                       9.1644428E-06, 7.5256819E-06, 7.0450071E-06, 7.0175175E-06, 6.7967953E-06, &
                                       5.9275452E-06, 4.4826864E-06, 3.1175764E-06, 2.4755370E-06, 2.7497674E-06, &
                                       2.3351376E-05, 2.5023929E-04, 4.0258401E-04, 3.9580811E-04, 3.5235096E-04, &
                                       3.4285950E-04, 4.2526209E-04]
      scav_coeff_mass(1:ncol, 10, 2) = [ &
                                       7.5335028E-05, 4.3935656E-05, 2.6325950E-05, 1.6456342E-05, 1.1086690E-05, &
                                       8.4099702E-06, 7.3245911E-06, 7.0370831E-06, 6.8749104E-06, 6.3065285E-06, &
                                       5.1618646E-06, 3.7982397E-06, 2.8388964E-06, 3.3371116E-06, 2.3043577E-05, &
                                       1.4169263E-04, 3.2300052E-04, 3.9310504E-04, 3.7321833E-04, 3.5162717E-04, &
                                       3.9570373E-04, 5.9692085E-04]
      scav_coeff_mass(1:ncol, 10, 3) = [ &
                                       5.5742505E-05, 3.2933162E-05, 2.0178400E-05, 1.3149277E-05, 9.4812690E-06, &
                                       7.7923218E-06, 7.1625847E-06, 6.8800873E-06, 6.3972453E-06, 5.4611997E-06, &
                                       4.2483097E-06, 3.3645392E-06, 5.6902336E-06, 2.9262298E-05, 1.1612305E-04, &
                                       2.5813040E-04, 3.5849826E-04, 3.7615875E-04, 3.6587389E-04, 4.0358739E-04, &
                                       5.8616078E-04, 1.1141649E-03]
      scav_coeff_mass(1:ncol, 10, 4) = [ &
                                       3.9661530E-05, 2.3954366E-05, 1.5273334E-05, 1.0641805E-05, 8.3659125E-06, &
                                       7.3761707E-06, 6.9063228E-06, 6.4002083E-06, 5.5773958E-06, 4.5557802E-06, &
                                       4.2988880E-06, 9.6367534E-06, 3.7078435E-05, 1.1111950E-04, 2.2518623E-04, &
                                       3.2357698E-04, 3.6565486E-04, 3.7655123E-04, 4.2652202E-04, 6.2461384E-04, &
                                       1.1915708E-03, 2.6566812E-03]
      scav_coeff_mass(1:ncol, 10, 5) = [ &
                                       2.7861705E-05, 1.7473000E-05, 1.1868066E-05, 9.0069450E-06, 7.6513762E-06, &
                                       6.9683530E-06, 6.3892151E-06, 5.6363022E-06, 4.9275818E-06, 5.9089737E-06, &
                                       1.4667872E-05, 4.5174097E-05, 1.1212815E-04, 2.0849267E-04, 2.9855833E-04, &
                                       3.5260152E-04, 3.8495088E-04, 4.5955770E-04, 6.9939975E-04, 1.3629470E-03, &
                                       3.0423535E-03, 6.9624439E-03]
      scav_coeff_mass(1:ncol, 11, 1) = [ &
                                       1.0227785E-04, 5.9328663E-05, 3.5227476E-05, 2.1668507E-05, 1.4229920E-05, &
                                       1.0494287E-05, 9.0492659E-06, 8.9212622E-06, 9.2527908E-06, 9.1913699E-06, &
                                       8.0974450E-06, 6.0598152E-06, 4.0237507E-06, 2.9422749E-06, 3.0846356E-06, &
                                       2.9637124E-05, 3.1285925E-04, 4.9697419E-04, 4.8839585E-04, 4.3249350E-04, &
                                       4.1345738E-04, 4.9750583E-04]
      scav_coeff_mass(1:ncol, 11, 2) = [ &
                                       8.1583280E-05, 4.7694013E-05, 2.8686394E-05, 1.8101602E-05, 1.2471760E-05, &
                                       9.8539437E-06, 9.0365904E-06, 9.0934876E-06, 9.1741654E-06, 8.5607184E-06, &
                                       7.0100663E-06, 5.0313034E-06, 3.5342087E-06, 3.9471213E-06, 2.8872922E-05, &
                                       1.7736938E-04, 4.0036510E-04, 4.8521558E-04, 4.5909934E-04, 4.2726902E-04, &
                                       4.6830140E-04, 6.8450663E-04]
      scav_coeff_mass(1:ncol, 11, 3) = [ &
                                       6.0430254E-05, 3.5825931E-05, 2.2112343E-05, 1.4663894E-05, 1.0942159E-05, &
                                       9.4396087E-06, 9.1064047E-06, 9.0699142E-06, 8.6140289E-06, 7.3953448E-06, &
                                       5.6682263E-06, 4.2908402E-06, 7.0097131E-06, 3.6555913E-05, 1.4490748E-04, &
                                       3.2035582E-04, 4.4296420E-04, 4.6276854E-04, 4.4529828E-04, 4.7948772E-04, &
                                       6.7466288E-04, 1.2529034E-03]
      scav_coeff_mass(1:ncol, 11, 4) = [ &
                                       4.3089621E-05, 2.6194871E-05, 1.6951267E-05, 1.2168127E-05, 1.0004993E-05, &
                                       9.2570923E-06, 9.0071510E-06, 8.5485126E-06, 7.5142733E-06, 6.0768631E-06, &
                                       5.5251192E-06, 1.2025173E-05, 4.6250276E-05, 1.3835645E-04, 2.7937103E-04, &
                                       3.9990233E-04, 4.4959849E-04, 4.5779641E-04, 5.0636283E-04, 7.1866754E-04, &
                                       1.3396452E-03, 2.9540811E-03]
      scav_coeff_mass(1:ncol, 11, 5) = [ &
                                       3.0420691E-05, 1.9332986E-05, 1.3485736E-05, 1.0668902E-05, 9.5000863E-06, &
                                       9.0007311E-06, 8.4670622E-06, 7.5480957E-06, 6.5415136E-06, 7.5771451E-06, &
                                       1.8351570E-05, 5.6285583E-05, 1.3940226E-04, 2.5845715E-04, 3.6876023E-04, &
                                       4.3293874E-04, 4.6669959E-04, 5.4360431E-04, 8.0234988E-04, 1.5299382E-03, &
                                       3.3808954E-03, 7.7081360E-03]
      scav_coeff_mass(1:ncol, 12, 1) = [ &
                                       1.1080475E-04, 6.4406221E-05, 3.8335775E-05, 2.3715529E-05, 1.5810478E-05, &
                                       1.2032807E-05, 1.0858323E-05, 1.1192844E-05, 1.1998009E-05, 1.2166364E-05, &
                                       1.0819968E-05, 8.0529470E-06, 5.1678713E-06, 3.5155814E-06, 3.4694541E-06, &
                                       3.6789490E-05, 3.8638738E-04, 6.1034052E-04, 6.0075583E-04, 5.2979458E-04, &
                                       4.9817835E-04, 5.8170076E-04]
      scav_coeff_mass(1:ncol, 12, 2) = [ &
                                       8.8437639E-05, 5.1820652E-05, 3.1289114E-05, 1.9937765E-05, 1.4053614E-05, &
                                       1.1549048E-05, 1.1092373E-05, 1.1603486E-05, 1.2015392E-05, 1.1374358E-05, &
                                       9.3352838E-06, 6.5879787E-06, 4.4029831E-06, 4.6661582E-06, 3.5604087E-05, &
                                       2.1920397E-04, 4.9256540E-04, 5.9639140E-04, 5.6327388E-04, 5.1850350E-04, &
                                       5.5396584E-04, 7.8379186E-04]
      scav_coeff_mass(1:ncol, 12, 3) = [ &
                                       6.5576317E-05, 3.9010936E-05, 2.4261634E-05, 1.6381040E-05, 1.2643950E-05, &
                                       1.1406169E-05, 1.1468976E-05, 1.1766317E-05, 1.1371214E-05, 9.8198254E-06, &
                                       7.4558361E-06, 5.4479161E-06, 8.5699421E-06, 4.5061811E-05, 1.7879487E-04, &
                                       3.9442179E-04, 5.4452192E-04, 5.6751933E-04, 5.4107906E-04, 5.6935040E-04, &
                                       7.7554919E-04, 1.4043225E-03]
      scav_coeff_mass(1:ncol, 12, 4) = [ &
                                       4.6862145E-05, 2.8679889E-05, 1.8845614E-05, 1.3937024E-05, 1.1953511E-05, &
                                       1.1536219E-05, 1.1587448E-05, 1.1213754E-05, 9.9353308E-06, 7.9852216E-06, &
                                       7.0419647E-06, 1.4851465E-05, 5.7018804E-05, 1.7054424E-04, 3.4392386E-04, &
                                       4.9152411E-04, 5.5084626E-04, 5.5553920E-04, 6.0075466E-04, 8.2579070E-04, &
                                       1.5011326E-03, 3.2691664E-03]
      scav_coeff_mass(1:ncol, 12, 5) = [ &
                                       3.3256178E-05, 2.1427563E-05, 1.5353780E-05, 1.2638205E-05, 1.1734722E-05, &
                                       1.1491981E-05, 1.1039721E-05, 9.9318832E-06, 8.5562546E-06, 9.6160927E-06, &
                                       2.2715685E-05, 6.9393061E-05, 1.7173298E-04, 3.1805351E-04, 4.5297385E-04, &
                                       5.2961593E-04, 5.6471015E-04, 6.4252685E-04, 9.1901892E-04, 1.7113447E-03, &
                                       3.7388646E-03, 8.4872283E-03]
      scav_coeff_mass(1:ncol, 13, 1) = [ &
                                       1.2010022E-04, 6.9939963E-05, 4.1724398E-05, 2.5951233E-05, 1.7546167E-05, &
                                       1.3739324E-05, 1.2888450E-05, 1.3769371E-05, 1.5141449E-05, 1.5603226E-05, &
                                       1.3991911E-05, 1.0392583E-05, 6.5159097E-06, 4.1843481E-06, 3.9023928E-06, &
                                       4.4620635E-05, 4.7134503E-04, 7.4543041E-04, 7.3631662E-04, 6.4739236E-04, &
                                       5.9953507E-04, 6.7967115E-04]
      scav_coeff_mass(1:ncol, 13, 2) = [ &
                                       9.5909453E-05, 5.6318797E-05, 3.4128908E-05, 2.1948447E-05, 1.5800047E-05, &
                                       1.3442019E-05, 1.3414107E-05, 1.4466222E-05, 1.5284491E-05, 1.4638140E-05, &
                                       1.2052446E-05, 8.4168828E-06, 5.4222745E-06, 5.4818249E-06, 4.3168068E-05, &
                                       2.6739198E-04, 6.0130423E-04, 7.2965705E-04, 6.8897893E-04, 6.2813027E-04, &
                                       6.5484390E-04, 8.9614674E-04]
      scav_coeff_mass(1:ncol, 13, 3) = [ &
                                       7.1186390E-05, 4.2485158E-05, 2.6612595E-05, 1.8272477E-05, 1.4539098E-05, &
                                       1.3621591E-05, 1.4157525E-05, 1.4861567E-05, 1.4560880E-05, 1.2644028E-05, &
                                       9.5493681E-06, 6.8019436E-06, 1.0348295E-05, 5.4765114E-05, 2.1805430E-04, &
                                       4.8159460E-04, 6.6567016E-04, 6.9347433E-04, 6.5609986E-04, 6.7549454E-04, &
                                       8.9037118E-04, 1.5688311E-03]
      scav_coeff_mass(1:ncol, 13, 4) = [ &
                                       5.0977372E-05, 3.1396961E-05, 2.0929887E-05, 1.5903813E-05, 1.4145397E-05, &
                                       1.4126432E-05, 1.4545286E-05, 1.4291699E-05, 1.2749285E-05, 1.0213734E-05, &
                                       8.8073902E-06, 1.8092276E-05, 6.9422136E-05, 2.0803447E-04, 4.1998165E-04, &
                                       6.0059365E-04, 6.7219993E-04, 6.7260307E-04, 7.1207086E-04, 9.4759711E-04, &
                                       1.6764319E-03, 3.5999548E-03]
      scav_coeff_mass(1:ncol, 13, 5) = [ &
                                       3.6356184E-05, 2.3730881E-05, 1.7428893E-05, 1.4851430E-05, 1.4272260E-05, &
                                       1.4345106E-05, 1.4007120E-05, 1.2697763E-05, 1.0902256E-05, 1.1980110E-05, &
                                       2.7747553E-05, 8.4590293E-05, 2.0954782E-04, 3.8838332E-04, 5.5314737E-04, &
                                       6.4520133E-04, 6.8167579E-04, 7.5866311E-04, 1.0509852E-03, 1.9073984E-03, &
                                       4.1137885E-03, 9.2914325E-03]
      scav_coeff_mass(1:ncol, 14, 1) = [ &
                                       1.3017396E-04, 7.5932455E-05, 4.5387641E-05, 2.8357187E-05, 1.9397740E-05, &
                                       1.5542620E-05, 1.5024650E-05, 1.6484636E-05, 1.8470214E-05, 1.9266448E-05, &
                                       1.7398017E-05, 1.2925152E-05, 7.9870425E-06, 4.9189243E-06, 4.3764819E-06, &
                                       5.2853756E-05, 5.6780803E-04, 9.0499195E-04, 8.9879648E-04, 7.8878561E-04, &
                                       7.2035162E-04, 7.9346229E-04]
      scav_coeff_mass(1:ncol, 14, 2) = [ &
                                       1.0400472E-04, 6.1186517E-05, 3.7193032E-05, 2.4103750E-05, 1.7654429E-05, &
                                       1.5438245E-05, 1.5860078E-05, 1.7492150E-05, 1.8758964E-05, 1.8129729E-05, &
                                       1.4979965E-05, 1.0401947E-05, 6.5362989E-06, 6.3666734E-06, 5.1429121E-05, &
                                       3.2186617E-04, 7.2810009E-04, 8.8817611E-04, 8.3976775E-04, 7.5927888E-04, &
                                       7.7334428E-04, 1.0230940E-03]
      scav_coeff_mass(1:ncol, 14, 3) = [ &
                                       7.7260577E-05, 4.6238237E-05, 2.9139122E-05, 2.0287939E-05, 1.6542923E-05, &
                                       1.5957880E-05, 1.6998539E-05, 1.8147357E-05, 1.7966393E-05, 1.5678316E-05, &
                                       1.1813329E-05, 8.2749557E-06, 1.2293126E-05, 6.5581381E-05, 2.6277878E-04, &
                                       5.8296778E-04, 8.0894174E-04, 8.4393679E-04, 7.9355541E-04, 8.0051376E-04, &
                                       1.0208557E-03, 1.7468410E-03]
      scav_coeff_mass(1:ncol, 14, 4) = [ &
                                       5.5425851E-05, 3.4321212E-05, 2.3156223E-05, 1.7988380E-05, 1.6460232E-05, &
                                       1.6864828E-05, 1.7684203E-05, 1.7574574E-05, 1.5767363E-05, 1.2617928E-05, &
                                       1.0726767E-05, 2.1681614E-05, 8.3424052E-05, 2.5104161E-04, 5.0854700E-04, &
                                       7.2926267E-04, 8.1663807E-04, 8.1208537E-04, 8.4295096E-04, 1.0858781E-03, &
                                       1.8659403E-03, 3.9440892E-03]
      scav_coeff_mass(1:ncol, 14, 5) = [ &
                                       3.9695768E-05, 2.6195476E-05, 1.9632774E-05, 1.7192427E-05, 1.6957125E-05, &
                                       1.7373266E-05, 1.7170528E-05, 1.5661225E-05, 1.3431421E-05, 1.4564229E-05, &
                                       3.3380726E-05, 1.0189344E-04, 2.5315989E-04, 4.7044478E-04, 6.7122527E-04, &
                                       7.8240138E-04, 8.2052540E-04, 8.9459557E-04, 1.1999898E-03, 2.1182990E-03, &
                                       4.5027397E-03, 1.0111084E-02]
      scav_coeff_mass(1:ncol, 15, 1) = [ &
                                       1.4102872E-04, 8.2381580E-05, 4.9315491E-05, 3.0908875E-05, 2.1315157E-05, &
                                       1.7351966E-05, 1.7118927E-05, 1.9121738E-05, 2.1702014E-05, 2.2838228E-05, &
                                       2.0742141E-05, 1.5434841E-05, 9.4649535E-06, 5.6752830E-06, 4.8807576E-06, &
                                       6.1135890E-05, 6.7519159E-04, 1.0915255E-03, 1.0920305E-03, 9.5772555E-04, &
                                       8.6369216E-04, 9.2529163E-04]
      scav_coeff_mass(1:ncol, 15, 2) = [ &
                                       1.1272353E-04, 6.6417260E-05, 4.0463339E-05, 2.6365015E-05, 1.9544514E-05, &
                                       1.7416292E-05, 1.8246201E-05, 2.0431469E-05, 2.2140862E-05, 2.1546222E-05, &
                                       1.7865830E-05, 1.2378838E-05, 7.6642543E-06, 7.2813044E-06, 6.0180477E-05, &
                                       3.8220265E-04, 8.7406430E-04, 1.0750446E-03, 1.0193695E-03, 9.1531856E-04, &
                                       9.1207752E-04, 1.1662723E-03]
      scav_coeff_mass(1:ncol, 15, 3) = [ &
                                       8.3793658E-05, 5.0254278E-05, 3.1806860E-05, 2.2362968E-05, 1.8546628E-05, &
                                       1.8248916E-05, 1.9763847E-05, 2.1345736E-05, 2.1294169E-05, 1.8661468E-05, &
                                       1.4057835E-05, 9.7556918E-06, 1.4327862E-05, 7.7344954E-05, 3.1281158E-04, &
                                       6.9929280E-04, 9.7670233E-04, 1.0222897E-03, 9.5684356E-04, 9.4720578E-04, &
                                       1.1688625E-03, 1.9387575E-03]
      scav_coeff_mass(1:ncol, 15, 4) = [ &
                                       6.0192095E-05, 3.7419340E-05, 2.5462840E-05, 2.0087495E-05, 1.8741543E-05, &
                                       1.9537216E-05, 2.0742429E-05, 2.0781608E-05, 1.8730715E-05, 1.4996658E-05, &
                                       1.2665522E-05, 2.5514295E-05, 9.8896628E-05, 2.9957414E-04, 6.1034770E-04, &
                                       8.7950301E-04, 9.8715241E-04, 9.7723734E-04, 9.9622103E-04, 1.2425555E-03, &
                                       2.0700430E-03, 4.2988945E-03]
      scav_coeff_mass(1:ncol, 15, 5) = [ &
                                       4.3241039E-05, 2.8759670E-05, 2.1864152E-05, 1.9509950E-05, 1.9584579E-05, &
                                       2.0327576E-05, 2.0261806E-05, 1.8569998E-05, 1.5938178E-05, 1.7217362E-05, &
                                       3.9496506E-05, 1.2121890E-04, 3.0270041E-04, 5.6500796E-04, 8.0898640E-04, &
                                       9.4390267E-04, 9.8429839E-04, 1.0530652E-03, 1.3678871E-03, 2.3442070E-03, &
                                       4.9024070E-03, 1.0935374E-02]
      scav_coeff_mass(1:ncol, 16, 1) = [ &
                                       1.5265715E-04, 8.9279445E-05, 5.3494687E-05, 3.3579958E-05, 2.3247062E-05, &
                                       1.9074724E-05, 1.9019262E-05, 2.1455763E-05, 2.4540929E-05, 2.5981443E-05, &
                                       2.3706130E-05, 1.7686673E-05, 1.0821211E-05, 6.4040610E-06, 5.4027420E-06, &
                                       6.9061124E-05, 7.9198513E-04, 1.3068683E-03, 1.3196213E-03, 1.1579575E-03, &
                                       1.0326699E-03, 1.0774028E-03]
      scav_coeff_mass(1:ncol, 16, 2) = [ &
                                       1.2205799E-04, 7.1999927E-05, 4.3919114E-05, 2.8691960E-05, 2.1396345E-05, &
                                       1.9251611E-05, 2.0382183E-05, 2.3023391E-05, 2.5115750E-05, 2.4564267E-05, &
                                       2.0437804E-05, 1.4167845E-05, 8.7166449E-06, 8.1809548E-06, 6.9144275E-05, &
                                       4.4750869E-04, 1.0395654E-03, 1.2929166E-03, 1.2313876E-03, 1.0996335E-03, &
                                       1.0736874E-03, 1.3273085E-03]
      scav_coeff_mass(1:ncol, 16, 3) = [ &
                                       9.0774592E-05, 5.4514054E-05, 3.4579392E-05, 2.4431235E-05, 2.0438478E-05, &
                                       2.0323041E-05, 2.2215585E-05, 2.4164073E-05, 2.4231896E-05, 2.1312630E-05, &
                                       1.6076525E-05, 1.1121661E-05, 1.6360749E-05, 8.9798706E-05, 3.6765215E-04, &
                                       8.3072783E-04, 1.1708061E-03, 1.2316764E-03, 1.1493153E-03, 1.1183847E-03, &
                                       1.3362431E-03, 2.1448788E-03]
      scav_coeff_mass(1:ncol, 16, 4) = [ &
                                       6.5256531E-05, 4.0655505E-05, 2.7785784E-05, 2.2094866E-05, 2.0827481E-05, &
                                       2.1920304E-05, 2.3444612E-05, 2.3614678E-05, 2.1361864E-05, 1.7132119E-05, &
                                       1.4474660E-05, 2.9456203E-05, 1.1560169E-04, 3.5334194E-04, 7.2562487E-04, &
                                       1.0528028E-03, 1.1864354E-03, 1.1712012E-03, 1.1746913E-03, 1.4195304E-03, &
                                       2.2890077E-03, 4.6613273E-03]
      scav_coeff_mass(1:ncol, 16, 5) = [ &
                                       4.6955055E-05, 3.1359233E-05, 2.4018432E-05, 2.1647297E-05, 2.1941032E-05, &
                                       2.2946276E-05, 2.2996531E-05, 2.1154117E-05, 1.8199956E-05, 1.9769087E-05, &
                                       4.5932412E-05, 1.4235732E-04, 3.5802417E-04, 6.7242304E-04, 9.6777116E-04, &
                                       1.1320760E-03, 1.1758800E-03, 1.2367601E-03, 1.5564856E-03, 2.5851350E-03, &
                                       5.3090489E-03, 1.1752415E-02]
      scav_coeff_mass(1:ncol, 17, 1) = [ &
                                       1.6503346E-04, 9.6608170E-05, 5.7907931E-05, 3.6345463E-05, 2.5149792E-05, &
                                       2.0634162E-05, 2.0599794E-05, 2.3299098E-05, 2.6737597E-05, 2.8408777E-05, &
                                       2.6015557E-05, 1.9475443E-05, 1.1942262E-05, 7.0610998E-06, 5.9312232E-06, &
                                       7.6203696E-05, 9.1547866E-04, 1.5515636E-03, 1.5843175E-03, 1.3927115E-03, &
                                       1.2300486E-03, 1.2517504E-03]
      scav_coeff_mass(1:ncol, 17, 2) = [ &
                                       1.3198590E-04, 7.7916436E-05, 4.7538433E-05, 3.1049054E-05, 2.3148062E-05, &
                                       2.0840982E-05, 2.2109651E-05, 2.5048990E-05, 2.7416947E-05, 2.6906627E-05, &
                                       2.2459852E-05, 1.5611409E-05, 9.6141523E-06, 9.0232889E-06, 7.7977696E-05, &
                                       5.1630837E-04, 1.2237707E-03, 1.5433919E-03, 1.4787393E-03, 1.3151704E-03, &
                                       1.2604942E-03, 1.5075312E-03]
      scav_coeff_mass(1:ncol, 17, 3) = [ &
                                       9.8182867E-05, 5.8995383E-05, 3.7423441E-05, 2.6436577E-05, 2.2125527E-05, &
                                       2.2037317E-05, 2.4154280E-05, 2.6354723E-05, 2.6512669E-05, 2.3389377E-05, &
                                       1.7689493E-05, 1.2264473E-05, 1.8296885E-05, 1.0258693E-04, 4.2634880E-04, &
                                       9.7649713E-04, 1.3920585E-03, 1.4744385E-03, 1.3737909E-03, 1.3164922E-03, &
                                       1.5245300E-03, 2.3651410E-03]
      scav_coeff_mass(1:ncol, 17, 4) = [ &
                                       7.0595256E-05, 4.3996019E-05, 3.0070260E-05, 2.3921647E-05, 2.2583017E-05, &
                                       2.3826770E-05, 2.5557899E-05, 2.5819447E-05, 2.3421952E-05, 1.8834929E-05, &
                                       1.6020017E-05, 3.3357219E-05, 1.3317178E-04, 4.1164219E-04, 8.5383998E-04, &
                                       1.2496985E-03, 1.4163463E-03, 1.3965176E-03, 1.3807469E-03, 1.6183531E-03, &
                                       2.5227140E-03, 5.0277459E-03]
      scav_coeff_mass(1:ncol, 17, 5) = [ &
                                       5.0802345E-05, 3.3938559E-05, 2.6007732E-05, 2.3473396E-05, 2.3847107E-05, &
                                       2.5007999E-05, 2.5132264E-05, 2.3181506E-05, 2.0022653E-05, 2.2060326E-05, &
                                       5.2492699E-05, 1.6494300E-04, 4.1858680E-04, 7.9234861E-04, 1.1480598E-03, &
                                       1.3484781E-03, 1.3975163E-03, 1.4478964E-03, 1.7672037E-03, 2.8406656E-03, &
                                       5.7182532E-03, 1.2549013E-02]
      scav_coeff_mass(1:ncol, 18, 1) = [ &
                                       1.7809283E-04, 1.0432826E-04, 6.2528146E-05, 3.9180876E-05, 2.6991825E-05, &
                                       2.1981027E-05, 2.1782000E-05, 2.4534397E-05, 2.8134010E-05, 2.9935111E-05, &
                                       2.7490369E-05, 2.0663967E-05, 1.2750716E-05, 7.6156230E-06, 6.4580016E-06, &
                                       8.2158532E-05, 1.0415600E-03, 1.8239835E-03, 1.8870028E-03, 1.6638010E-03, &
                                       1.4575097E-03, 1.4494069E-03]
      scav_coeff_mass(1:ncol, 18, 2) = [ &
                                       1.4245458E-04, 8.4133038E-05, 5.1294933E-05, 3.3407479E-05, 2.4758195E-05, &
                                       2.2119248E-05, 2.3329529E-05, 2.6370169E-05, 2.8873834E-05, 2.8393021E-05, &
                                       2.3775880E-05, 1.6603543E-05, 1.0302431E-05, 9.7746819E-06, 8.6287555E-05, &
                                       5.8645576E-04, 1.4240812E-03, 1.8260933E-03, 1.7627120E-03, 1.5636257E-03, &
                                       1.4738323E-03, 1.7074241E-03]
      scav_coeff_mass(1:ncol, 18, 3) = [ &
                                       1.0597718E-04, 6.3668168E-05, 4.0309408E-05, 2.8339685E-05, 2.3548151E-05, &
                                       2.3301897E-05, 2.5454118E-05, 2.7759672E-05, 2.7963748E-05, 2.4732441E-05, &
                                       1.8776725E-05, 1.3109722E-05, 2.0048766E-05, 1.1525472E-04, 4.8739560E-04, &
                                       1.1344722E-03, 1.6394420E-03, 1.7512058E-03, 1.6317148E-03, 1.5428874E-03, &
                                       1.7343491E-03, 2.5986092E-03]
      scav_coeff_mass(1:ncol, 18, 4) = [ &
                                       7.6173633E-05, 4.7408949E-05, 3.2276505E-05, 2.5509902E-05, 2.3922694E-05, &
                                       2.5138196E-05, 2.6933770E-05, 2.7231686E-05, 2.4754665E-05, 1.9979094E-05, &
                                       1.7205327E-05, 3.7062502E-05, 1.5109374E-04, 4.7323486E-04, 9.9330346E-04, &
                                       1.4691059E-03, 1.6770700E-03, 1.6542876E-03, 1.6156121E-03, 1.8396035E-03, &
                                       2.7701166E-03, 5.3933893E-03]
      scav_coeff_mass(1:ncol, 18, 5) = [ &
                                       5.4747810E-05, 3.6456097E-05, 2.7773864E-05, 2.4904672E-05, 2.5189024E-05, &
                                       2.6371427E-05, 2.6512596E-05, 2.4500495E-05, 2.1275831E-05, 2.3967519E-05, &
                                       5.8956393E-05, 1.8842082E-04, 4.8329990E-04, 9.2339659E-04, 1.3488678E-03, &
                                       1.5930743E-03, 1.6500033E-03, 1.6874760E-03, 2.0004287E-03, 3.1093862E-03, &
                                       6.1243845E-03, 1.3310020E-02]
      scav_coeff_mass(1:ncol, 19, 1) = [ &
                                       1.9168898E-04, 1.1235384E-04, 6.7304105E-05, 4.2054330E-05, 2.8750759E-05, &
                                       2.3094934E-05, 2.2540989E-05, 2.5126532E-05, 2.8681312E-05, 3.0499569E-05, &
                                       2.8067120E-05, 2.1200355E-05, 1.3214927E-05, 8.0533770E-06, 6.9774789E-06, &
                                       8.6583777E-05, 1.1646868E-03, 2.1192455E-03, 2.2252210E-03, 1.9702127E-03, &
                                       1.7144678E-03, 1.6695920E-03]
      scav_coeff_mass(1:ncol, 19, 2) = [ &
                                       1.5334703E-04, 9.0580621E-05, 5.5146679E-05, 3.5739753E-05, 2.6204950E-05, &
                                       2.3063287E-05, 2.4011312E-05, 2.6944617E-05, 2.9431675E-05, 2.8961988E-05, &
                                       2.4329150E-05, 1.7103055E-05, 1.0758473E-05, 1.0412955E-05, 9.3654628E-05, &
                                       6.5511498E-04, 1.6355263E-03, 2.1374149E-03, 2.0815400E-03, 1.8441574E-03, &
                                       1.7129744E-03, 1.9257231E-03]
      scav_coeff_mass(1:ncol, 19, 3) = [ &
                                       1.1407056E-04, 6.8480203E-05, 4.3204060E-05, 3.0115962E-05, 2.4682692E-05, &
                                       2.4087777E-05, 2.6076111E-05, 2.8328551E-05, 2.8527190E-05, 2.5285278E-05, &
                                       1.9292986E-05, 1.3625821E-05, 2.1542790E-05, 1.2725796E-04, 5.4866089E-04, &
                                       1.3007199E-03, 1.9091089E-03, 2.0595752E-03, 1.9218465E-03, 1.7967050E-03, &
                                       1.9644533E-03, 2.8426160E-03]
      scav_coeff_mass(1:ncol, 19, 4) = [ &
                                       8.1929155E-05, 5.0855074E-05, 3.4376542E-05, 2.6834398E-05, 2.4817482E-05, &
                                       2.5817110E-05, 2.7524699E-05, 2.7796659E-05, 2.5305242E-05, 2.0517446E-05, &
                                       1.7982840E-05, 4.0419783E-05, 1.6870077E-04, 5.3622996E-04, 1.1407564E-03, &
                                       1.7074572E-03, 1.9659209E-03, 1.9428990E-03, 1.8781831E-03, 2.0818789E-03, &
                                       3.0283368E-03, 5.7514477E-03]
      scav_coeff_mass(1:ncol, 19, 5) = [ &
                                       5.8746038E-05, 3.8880166E-05, 2.9289491E-05, 2.5911305E-05, 2.5929870E-05, &
                                       2.6990919E-05, 2.7085421E-05, 2.5058062E-05, 2.1908225E-05, 2.5413902E-05, &
                                       6.5081051E-05, 2.1201730E-04, 5.5038246E-04, 1.0627157E-03, 1.5669596E-03, &
                                       1.8631447E-03, 1.9314717E-03, 1.9541224E-03, 2.2544742E-03, 3.3879365E-03, &
                                       6.5195931E-03, 1.4017044E-02]
      scav_coeff_mass(1:ncol, 20, 1) = [ &
                                       2.0552655E-04, 1.2051263E-04, 7.2136228E-05, 4.4911228E-05, 3.0402518E-05, &
                                       2.3975475E-05, 2.2896877E-05, 2.5113453E-05, 2.8430207E-05, 3.0156386E-05, &
                                       2.7791574E-05, 2.1113051E-05, 1.3345911E-05, 8.3741706E-06, 7.4837191E-06, &
                                       8.9237650E-05, 1.2781326E-03, 2.4281236E-03, 2.5913814E-03, 2.3062741E-03, &
                                       1.9964904E-03, 1.9083615E-03]
      scav_coeff_mass(1:ncol, 20, 2) = [ &
                                       1.6442758E-04, 9.7122332E-05, 5.9016209E-05, 3.8006511E-05, 2.7476266E-05, &
                                       2.3683194E-05, 2.4184145E-05, 2.6816457E-05, 2.9142336E-05, 2.8662710E-05, &
                                       2.4156148E-05, 1.7129508E-05, 1.0987776E-05, 1.0925914E-05, 9.9666822E-05, &
                                       7.1884859E-04, 1.8502739E-03, 2.4691078E-03, 2.4286202E-03, 2.1516940E-03, &
                                       1.9736889E-03, 2.1581908E-03]
      scav_coeff_mass(1:ncol, 20, 3) = [ &
                                       1.2228976E-04, 7.3332540E-05, 4.6054712E-05, 3.1744245E-05, 2.5532097E-05, &
                                       2.4417831E-05, 2.6058928E-05, 2.8109467E-05, 2.8251491E-05, 2.5087365E-05, &
                                       1.9263054E-05, 1.3820907E-05, 2.2721135E-05, 1.3798822E-04, 6.0738545E-04, &
                                       1.4691289E-03, 2.1932944E-03, 2.3925075E-03, 2.2385740E-03, 2.0733403E-03, &
                                       2.2104176E-03, 3.0915805E-03]
      scav_coeff_mass(1:ncol, 20, 4) = [ &
                                       8.7742108E-05, 5.4269397E-05, 3.6341435E-05, 2.7892516E-05, 2.5285579E-05, &
                                       2.5897873E-05, 2.7375137E-05, 2.7560793E-05, 2.5113530E-05, 2.0476512E-05, &
                                       1.8350638E-05, 4.3283073E-05, 1.8518179E-04, 5.9802566E-04, 1.2909976E-03, &
                                       1.9577753E-03, 2.2759211E-03, 2.2564148E-03, 2.1634953E-03, 2.3404343E-03, &
                                       3.2914209E-03, 6.0917602E-03]
      scav_coeff_mass(1:ncol, 20, 5) = [ &
                                       6.2720054E-05, 4.1174646E-05, 3.0547220E-05, 2.6507692E-05, 2.6100459E-05, &
                                       2.6907617E-05, 2.6894701E-05, 2.4892848E-05, 2.1942605E-05, 2.6367369E-05, &
                                       7.0604466E-05, 2.3472897E-04, 6.1724958E-04, 1.2055970E-03, 1.7959987E-03, &
                                       2.1519940E-03, 2.2358672E-03, 2.2425669E-03, 2.5241881E-03, 3.6697062E-03, &
                                       6.8924256E-03, 1.4646616E-02]
      scav_coeff_mass(1:ncol, 21, 1) = [ &
                                       2.1908395E-04, 1.2850022E-04, 7.6848478E-05, 4.7655422E-05, 3.1906320E-05, &
                                       2.4626989E-05, 2.2896333E-05, 2.4582640E-05, 2.7502284E-05, 2.9043693E-05, &
                                       2.6790039E-05, 2.0489711E-05, 1.3185023E-05, 8.5853992E-06, 7.9658604E-06, &
                                       9.0002160E-05, 1.3745584E-03, 2.7363928E-03, 2.9711936E-03, 2.6599098E-03, &
                                       2.2937466E-03, 2.1573088E-03]
      scav_coeff_mass(1:ncol, 21, 2) = [ &
                                       1.7528008E-04, 1.0351635E-04, 6.2766829E-05, 4.0139373E-05, 2.8554465E-05, &
                                       2.4005215E-05, 2.3915669E-05, 2.6090122E-05, 2.8134534E-05, 2.7625441E-05, &
                                       2.3362038E-05, 1.6746624E-05, 1.1014864E-05, 1.1306745E-05, 1.0395807E-04, &
                                       7.7384755E-04, 2.0575174E-03, 2.8071875E-03, 2.7908872E-03, 2.4753034E-03, &
                                       2.2468198E-03, 2.3964051E-03]
      scav_coeff_mass(1:ncol, 21, 3) = [ &
                                       1.3032880E-04, 7.8050572E-05, 4.8769346E-05, 3.3190354E-05, 2.6108970E-05, &
                                       2.4346696E-05, 2.5494384E-05, 2.7220826E-05, 2.7262641E-05, 2.4248640E-05, &
                                       1.8762751E-05, 1.3731080E-05, 2.3540348E-05, 1.4681551E-04, 6.6029641E-04, &
                                       1.6313132E-03, 2.4795467E-03, 2.7369353E-03, 2.5703267E-03, 2.3629731E-03, &
                                       2.4633539E-03, 3.3358486E-03]
      scav_coeff_mass(1:ncol, 21, 4) = [ &
                                       9.3401328E-05, 5.7538256E-05, 3.8123289E-05, 2.8686795E-05, 2.5372597E-05, &
                                       2.5463114E-05, 2.6594594E-05, 2.6643724E-05, 2.4288488E-05, 1.9936544E-05, &
                                       1.8340396E-05, 4.5515510E-05, 1.9961925E-04, 6.5535631E-04, 1.4367372E-03, &
                                       2.2090285E-03, 2.5946042E-03, 2.5830911E-03, 2.4612501E-03, 2.6058505E-03, &
                                       3.5491240E-03, 6.3995666E-03]
      scav_coeff_mass(1:ncol, 21, 5) = [ &
                                       6.6535252E-05, 4.3279183E-05, 3.1541287E-05, 2.6732520E-05, 2.5776293E-05, &
                                       2.6223505E-05, 2.6053515E-05, 2.4110183E-05, 2.1455832E-05, 2.6829258E-05, &
                                       7.5249797E-05, 2.5534814E-04, 6.8050955E-04, 1.3452775E-03, 2.0259457E-03, &
                                       2.4478821E-03, 2.5515784E-03, 2.5422126E-03, 2.7996007E-03, 3.9435675E-03, &
                                       7.2265035E-03, 1.5168506E-02]
      scav_coeff_mass(1:ncol, 22, 1) = [ &
                                       2.3156777E-04, 1.3585274E-04, 8.1170971E-05, 5.0136448E-05, 3.3192574E-05, &
                                       2.5043735E-05, 2.2592101E-05, 2.3642986E-05, 2.6054297E-05, 2.7343561E-05, &
                                       2.5231953E-05, 1.9449422E-05, 1.2788074E-05, 8.6945987E-06, 8.4040313E-06, &
                                       8.8889600E-05, 1.4468730E-03, 3.0251956E-03, 3.3432472E-03, 3.0119083E-03, &
                                       2.5902814E-03, 2.4029527E-03]
      scav_coeff_mass(1:ncol, 22, 2) = [ &
                                       1.8527095E-04, 1.0939330E-04, 6.6187340E-05, 4.2028279E-05, 2.9402523E-05, &
                                       2.4053824E-05, 2.3287504E-05, 2.4898399E-05, 2.6576247E-05, 2.6023309E-05, &
                                       2.2088528E-05, 1.6040646E-05, 1.0871527E-05, 1.1548706E-05, 1.0624810E-04, &
                                       8.1631053E-04, 2.2440446E-03, 3.1318983E-03, 3.1482536E-03, 2.7974770E-03, &
                                       2.5176235E-03, 2.6272062E-03]
      scav_coeff_mass(1:ncol, 22, 3) = [ &
                                       1.3772091E-04, 8.2365801E-05, 5.1202591E-05, 3.4393205E-05, 2.6418529E-05, &
                                       2.3938065E-05, 2.4497897E-05, 2.5816067E-05, 2.5727036E-05, 2.2916272E-05, &
                                       1.7894237E-05, 1.3405421E-05, 2.3969208E-05, 1.5314899E-04, 7.0387457E-04, &
                                       1.7770311E-03, 2.7508573E-03, 3.0734041E-03, 2.8989422E-03, 2.6499019E-03, &
                                       2.7093267E-03, 3.5612203E-03]
      scav_coeff_mass(1:ncol, 22, 4) = [ &
                                       9.8582809E-05, 6.0483578E-05, 3.9640602E-05, 2.9207916E-05, 2.5129638E-05, &
                                       2.4615656E-05, 2.5324222E-05, 2.5202721E-05, 2.2975258E-05, 1.9005615E-05, &
                                       1.8001614E-05, 4.6994912E-05, 2.1106425E-04, 7.0451213E-04, 1.5689061E-03, &
                                       2.4462674E-03, 2.9037811E-03, 2.9048456E-03, 2.7551834E-03, 2.8634449E-03, &
                                       3.7864224E-03, 6.6551643E-03]
      scav_coeff_mass(1:ncol, 22, 5) = [ &
                                       6.9981710E-05, 4.5093472E-05, 3.2250139E-05, 2.6627101E-05, 2.5050415E-05, &
                                       2.5069424E-05, 2.4709959E-05, 2.2849985E-05, 2.0553105E-05, 2.6820626E-05, &
                                       7.8740479E-05, 2.7254654E-04, 7.3615275E-04, 1.4731772E-03, 2.2431651E-03, &
                                       2.7338469E-03, 2.8609875E-03, 2.8365528E-03, 3.0653482E-03, 4.1933884E-03, &
                                       7.5001904E-03, 1.5545671E-02]

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_impc_scav_dust_init

! ----------------------------------------------------------------------
   SUBROUTINE ukca_impc_scav_dust_dealloc()
! ----------------------------------------------------------------------
! Purpose:
! -------
! Deallocate scavenging coefficient lookup tables
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UKCA_IMPC_SCAV_DUST_DEALLOC'

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

      IF (ALLOCATED(scav_coeff_num)) DEALLOCATE (scav_coeff_num)
      IF (ALLOCATED(scav_coeff_mass)) DEALLOCATE (scav_coeff_mass)

      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
      RETURN
   END SUBROUTINE ukca_impc_scav_dust_dealloc

END MODULE ukca_impc_scav_dust_mod
