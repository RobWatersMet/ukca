! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: BOX_MODEL                     -------------------------
!
!    Purpose: Output chemical fluxes and rates
!             Currently done in a method not via UKCA API - needs fixing!  
!
!    -------------------------------------------------------------------
MODULE box_output_chem_diags_mod
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: box_output_chem_diags

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BOX_OUTPUT_CHEM_DIAGS_MOD'
  
CONTAINS

  SUBROUTINE box_output_chem_diags(timestep, flux_unit, rate_unit)

    ! Take flux and rate arrays directly from ASAD_MOD for the moment
    ! Highly bad thing to do - really need to have these passed through UKCA API

    ! prk is in units of molecules.cm^-3.s^-1
    ! rk is in units of s^-1
    ! speci is the list of species in ASAD
    ! nspi is the ordering of species within the *rk arrays, of size (jpnr,jpmsp[==6])
    
    USE asad_mod, ONLY: rk, prk, jpnr, speci, nspi, jpmsp, jpspec
    
    USE ukca_api_mod, ONLY: ukca_maxlen_fieldname

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: timestep
    INTEGER, INTENT(IN) :: flux_unit
    INTEGER, INTENT(IN) :: rate_unit
    
    CHARACTER(LEN=*), PARAMETER :: RoutineName='BOX_OUTPUT_CHEM_DIAGS'

    ! Format string for outputting of box model
    CHARACTER(LEN=80) :: out_format, head_format, col_format
    CHARACTER(LEN=15000) :: header

    ! local speci array to cope with 0 locations
    CHARACTER(LEN=10) :: zero_speci(0:jpspec)

    ! string to print describing reactions
    CHARACTER(LEN=8), PARAMETER :: prodreact(jpmsp)=&
         (/'  PROD 1', '  PROD 2', ' REACT 1',' REACT 2',' REACT 3',' REACT 4'/) 
    
    ! number of gridpoints = will be 1 for box model
    INTEGER, PARAMETER :: n_points=1

    ! only do header on first pass through
    LOGICAL, SAVE :: l_first=.TRUE.

    INTEGER :: i,j

    ! output format for the fluxes/rates
    WRITE(out_format,  '(A5, I0, A18)')  '(I10,',         jpnr, '(",",1X,ES15.6E3))'
    WRITE(head_format, '(A9, I0, A13)')  '("# ",A8,',     jpnr, '(",",1X,A15))'
    WRITE(col_format,  '(A13, I0, A13)') '("# COL ",I4,', jpnr, '(",",1X,I15))'
    
    IF (l_first) THEN
       ! first line of header
       WRITE(flux_unit, '(A)')                                                  &
            '# Timestep, UKCA reaction fluxes (units = molecules.cm^-3.s^-1)'
       WRITE(rate_unit, '(A)') '# Timestep, UKCA reaction rates (units = s^-1)'

       ! set-up species array to write out
       zero_speci(0) = '         '
       zero_speci(1:jpspec) = speci(:)

       ! write column numbers
       WRITE(header, TRIM(ADJUSTL(col_format))) 1,(i+1, i=1,jpnr)
       WRITE(flux_unit, '(A)') TRIM(ADJUSTL(header))
       WRITE(rate_unit, '(A)') TRIM(ADJUSTL(header))
       
       ! now write full header information for each file (identical)
       DO j=1,jpmsp
          WRITE(header, TRIM(ADJUSTL(head_format))) TRIM(ADJUSTL(prodreact(j))),&
               (TRIM(ADJUSTL(zero_speci(nspi(i,j)))), i=1,jpnr)
          WRITE(flux_unit, '(A)') TRIM(ADJUSTL(header))
          WRITE(rate_unit, '(A)') TRIM(ADJUSTL(header))
       END DO
       
       l_first=.FALSE.
    END IF

    ! write out fluxes from prk array
    WRITE(flux_unit, TRIM(ADJUSTL(out_format))) timestep, prk(1,:)

    ! write out rates from rk array
    WRITE(rate_unit, TRIM(ADJUSTL(out_format))) timestep, rk(1,:)
       
    RETURN
  END SUBROUTINE box_output_chem_diags

END MODULE box_output_chem_diags_mod
