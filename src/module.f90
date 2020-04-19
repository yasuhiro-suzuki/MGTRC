!=module.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!==Reference
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!

MODULE kind_spec
  INTEGER, PARAMETER :: DP =  SELECTED_REAL_KIND(15)
END MODULE kind_spec

MODULE param1
  USE kind_spec
  IMPLICIT NONE
  REAL(DP), PARAMETER :: pi  =  3.141592653589793238462643383279502884197_DP, &
    &                    pi2 =  pi  + pi,                                     &
    &                    pi4 =  pi2 + pi2,                                    &
    &                    mu0 =  pi4 * 1.0E-07_DP,                             &
    &                    c12 =  0.5_DP,                                       &
    &                    c13 =  1.0_DP / 3.0_DP,                              &
    &                    c14 =  0.25_DP
  CHARACTER(LEN=5) :: ver_info =  ''
END MODULE param1

MODULE fline_mod
  USE kind_spec
  IMPLICIT NONE
  CHARACTER(LEN=20) :: mode      = '2D', &
    &                  def_start = ''
  LOGICAL :: lupdown = .false., &
    &        lflxqnt = .false., &
    &        lpoin   = .false., &
    &        ltext   = .false., &
    &        lfigout = .false., &
    &        lcontb
  INTEGER :: mr     =  30,  &
    &        nstep  =  32,  &
    &        mcirc  =  100, &
    &        ntheta =  150
  REAL(DP) :: h_in        =  0.01_DP,    &
    &         lc_in       =  1.0e+03_DP, &
    &         drflx       =  0.02_DP,    &
    &         dzflx       =  0.0_DP,     &
    &         dpflx       =  0.0_DP,     &
    &         rstart      =  0.0_DP,     &
    &         zstart      =  0.0_DP,     &
    &         pstart      =  0.0_DP,     &
    &         sstart      =  0.0_DP,     &
    &         sout        =  1.0_DP,     &
    &         rout        =  0.0_DP,     &
    &         zout        =  0.0_DP,     &
    &         pcros_in(8) =  0.0_DP
END MODULE fline_mod

MODULE axis_mod
  USE kind_spec
  IMPLICIT NONE
  LOGICAL :: lrmagax
  REAL(DP) :: rax0    =  0.0_DP, &
    &         zax0    =  0.0_DP, &
    &         raxis,             &
    &         zaxis,             &
    &         bpaxis,            &
    &         baxis,             &
    &         paxis,             &
    &         rax_av,            &
    &         zax_av,            &
    &         bax_av,            &
    &         bpax_av,           &
    &         rmagax,            &
    &         fbvadj,            &
    &         btor0,             &
    &         rbtor0
END MODULE axis_mod

MODULE file_mod
  IMPLICIT NONE
  LOGICAL :: lfile
  CHARACTER(LEN=300) :: plot_file_name = 'default', &
#ifdef POINXY
    &                   xy_plot_file   = '',        &
#endif
#ifdef LCOUT
    &                   lc_plot_file   = '',        &
#endif
    &                   prof_file      = ''
END MODULE file_mod
