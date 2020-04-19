!=file_close.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! File close of input and output files
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
! None 
!
!==TODO
!
! None
!

SUBROUTINE file_close

  USE fline_mod,             ONLY : ltext
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
    &                               flx_file

  IMPLICIT NONE


  CLOSE(25)

  IF(lflux) CLOSE(26)

  IF(ltext) CLOSE(60)

#ifdef POINXY
  CLOSE(51)
#endif
#ifdef LCOUT
  CLOSE(52)
#endif


END SUBROUTINE file_close
