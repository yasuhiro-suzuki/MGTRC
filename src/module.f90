!> @file module.f90
!------------------------------------------------------------------------------
!
! MODULE: kind_spec
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> @brief
!> Prescribe precision.
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE kind_spec
  INTEGER, PARAMETER :: DP =  SELECTED_REAL_KIND(15) !< Double Precision
END MODULE kind_spec
!------------------------------------------------------------------------------
!
! MODULE: param1
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> Definition of constants
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE param1
  USE kind_spec
  IMPLICIT NONE
  REAL(DP), PARAMETER :: pi  =  3.141592653589793238462643383279502884197_DP, & !< Pi
    &                    pi2 =  pi  + pi,                                     & !< 2 * Pi
    &                    pi4 =  pi2 + pi2,                                    & !< 4 * Pi
    &                    mu0 =  pi4 * 1.0E-07_DP,                             & !< Magnetic permeability for the vacuum
    &                    c12 =  0.5_DP,                                       & !< 1/2
    &                    c13 =  1.0_DP / 3.0_DP,                              & !< 1/3
    &                    c14 =  0.25_DP                                         !< 1/4
  CHARACTER(LEN=5) :: ver_info =  '1.0.0'                                       !< Version of code
END MODULE param1
!------------------------------------------------------------------------------
!
! MODULE: fline_mod
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> Prescribe control parameters for the field line tracer.
!
!> @details
!> Poloidal cross section to start field line tracong is defined by
!! \f[
!! \phi = \frac{2 \pi}{M} \times \mathrm{pstart}
!! \f]
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE fline_mod
  USE kind_spec
  IMPLICIT NONE
  CHARACTER(LEN=20) :: mode      = '2D', & !< Switch of the integral variable
    &                  def_start = ''      !< Definition of starting points for the filed line trace
  LOGICAL :: lupdown = .false., &          !< Logical switch of the direction of the integral
    &        lpoin   = .false., &          !< Logical switch to calculate Poincare plots
    &        lflxqnt = .false., &          !< Logical switch to calculate flux surface quantities
    &        lfigout = .false., &          !< Logical switch to output the graphic file
    &        lcontb  = .false.             !< Logical switch to draw the magnetic field strength
  INTEGER :: mr     =  30,  &              !< Number of flux surfaces in the field line tracing 
    &        nstep  =  32,  &              !< Number of steps in the one toroidal field period (mode: 2D)
    &        mcirc  =  100, &              !< Number of toroidal circuits (mode: 2D)
    &        ntheta =  150                 !< Number of points of Poincare plots in the subroutine, mgvmec
  REAL(DP) :: h_in        =  0.01_DP,    & !< Step size [m] for the field line tracer (mode: 3D)
    &         lc_in       =  1.0e+03_DP, & !< Length [m] of the integral (mode: 3D)
    &         drflx       =  0.02_DP,    & !< Step size [m] to trace the field line along R dicrection
    &         dzflx       =  0.0_DP,     & !< Step size [m] to trace the field line along Z dicrection
    &         dpflx       =  0.0_DP,     & !< Step size [rad] to trace the field line along phi dicrection
    &         rstart      =  0.0_DP,     & !< R [m] to start the field line tracer
    &         zstart      =  0.0_DP,     & !< Z [m] to start the field line tracer
    &         pstart      =  0.0_DP,     & !< Normalized phi to start the field line tracer
    &         sstart      =  0.0_DP,     & !< Normalized flux to start the field line tracer
    &         rout        =  0.0_DP,     & !< R [m] to end the field line tracer
    &         zout        =  0.0_DP,     & !< Z [m] to end the field line tracer
    &         sout        =  0.0_DP,     & !< Normalized flux to end the field line tracer
    &         pcros_in(8) =  0.0_DP        !< Additional toroidal angle to calculate Poincare plots up to 8
END MODULE fline_mod
!------------------------------------------------------------------------------
!
! MODULE: axis_mod
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> Prescribe to search the magnetic axis, and store axis positions, field and so on
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE axis_mod
  USE kind_spec
  IMPLICIT NONE
  LOGICAL :: lrmagax               !<
  REAL(DP) :: rax0    =  0.0_DP, & !< R [m] of the initial guess
    &         zax0    =  0.0_DP, & !< Z [m] of the initial guess
    &         raxis,             & !< R [m] of the axis
    &         zaxis,             & !< Z [m] of the axis
    &         bpaxis,            & !< Toroidal magnetic field [T] on the axis
    &         baxis,             & !< Magnetic field [T] on the axis
    &         rax_av,            & !< \<R\> [m] of the toroidally averaged axis
    &         zax_av,            & !< \<Z\> [m] of the toroidally averaged axis
    &         bpax_av,           & !< Toroidally averaged \<B_p\> [T]
    &         bax_av,            & !< Toroidally averaged \<B\> [T]
    &         rmagax,            & !<
    &         fbvadj,            & !<
    &         btor0,             & !<
    &         rbtor0               !<
END MODULE axis_mod
!------------------------------------------------------------------------------
!
! MODULE: file_mod
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> Define names of output files
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
MODULE file_mod
  IMPLICIT NONE
  LOGICAL :: lfile                                    !< Logical switch to use OPEN/CLOSE statements
  CHARACTER(LEN=300) :: plot_file_name = 'default', & !< Base name of fiels for Poincare plots
#ifdef POINXY
    &                   xy_plot_file   = '',        & !< File name of Poincare plots on X-Y plane
#endif
#ifdef LCOUT
    &                   lc_plot_file   = '',        & !< File name of the connection lenght plot
#endif
    &                   prof_file      = ''           !< File name of radial profiles of the flux surface quantity
END MODULE file_mod
