!> @file file_open.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: file_open
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> @brief
!>
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
!=file_open.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! File open of input and output files
!
! Device no. 25 : magnetic field
! Device no. 26 : toroidal flux
! Device no. 51 : Poincare plot on X-Y plane (Option)
! Device no. 52 : L_C on a line (Option)
! Device no. 60 : radial profile of flux surface quantities
!
!==Reference
!
! None
!
!==Error Handlings
!
! None
!
!==Known Bugs
!
! None
!
!==Note
!
! For using XL Fortran, the user should take care the compatibility
! of open statement 
!
!==TODO
!
! Implementation of XSD, netCDF and HDF5 file systems
!

SUBROUTINE file_open

  USE fline_mod,             ONLY : lflxqnt
  USE file_mod,              ONLY : plot_file_name, &
#ifdef POINXY
    &                               xy_plot_file,   &
#endif
#ifdef LCOUT
    &                               lc_plot_file,   &
#endif
    &                               prof_file
  USE cylindrical_coord_mod, ONLY : mag_file,       &
    &                               lflux,          &
    &                               flx_file,       &
    &                               flx_form
  USE vessel_mod,            ONLY : lvessel,        &
    &                               vessel_file

  IMPLICIT NONE


  OPEN(25, FILE=mag_file,     FORM='unformatted', STATUS='old')

  IF(lflux)THEN
    SELECT CASE(TRIM(flx_form))
      CASE('xss')
        OPEN(26, FILE=flx_file, FORM='unformatted', STATUS='old')
      CASE('flx')
        OPEN(26, FILE=flx_file, FORM='unformatted', STATUS='old') 
      CASE('eqdsk')
        OPEN(26, FILE=flx_file, FORM='formatted', STATUS='old') 
      CASE('eqdata')
        OPEN(26, FILE=flx_file, FORM='formatted', STATUS='old') 
      CASE DEFAULT
        OPEN(26, FILE=flx_file, FORM='unformatted', STATUS='old')
    END SELECT
  END IF

  IF(lvessel)THEN
    OPEN(30, FILE=vessel_file, FORM='formatted', STATUS='unknown')
  END IF
 
  IF(lflxqnt)THEN
    IF(prof_file == '') prof_file =  TRIM(plot_file_name) // '.prof'
    OPEN(60, FILE=prof_file, FORM='formatted', STATUS='unknown')
  END IF

#ifdef POINXY
  IF(xy_plot_file == '') xy_plot_file =  plot_file_name // '.xy'
  OPEN(51, FILE=xy_plot_file, FORM='formatted', STATUS='unknown')
#endif
#ifdef LCOUT
  IF(lc_plot_file == '') lc_plot_file =  TRIM(plot_file_name) // '.lc'
  OPEN(52, FILE=lc_plot_file, FORM='formatted', STATUS='unknown')
#endif


END SUBROUTINE file_open
