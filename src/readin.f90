!> @file readin.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: readin
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
SUBROUTINE readin

  USE kind_spec
  USE file_mod,              ONLY : lfile,          &
    &                               plot_file_name, &
#ifdef POINXY
    &                               xy_plot_file,   &
#endif
#ifdef LCOUT
    &                               lc_plot_file,   &
#endif
    &                               prof_file
  USE axis_mod,              ONLY : lrmagax,        &
    &                               rax0,           &
    &                               zax0,           &
    &                               rmagax,         &
    &                               fbvadj,         &
    &                               btor0,          &
    &                               rbtor0
  USE fline_mod,             ONLY : mode,           &
    &                               def_start,      &
    &                               lupdown,        &
    &                               lflxqnt,        &
    &                               lpoin,          &
    &                               lfigout,        &
    &                               lcontb,         &
    &                               mr,             &
    &                               nstep,          &
    &                               mcirc,          &
    &                               ntheta,         &
    &                               h_in,           &
    &                               lc_in,          &
    &                               drflx,          &
    &                               dzflx,          &
    &                               dpflx,          &
    &                               rstart,         &
    &                               zstart,         &
    &                               pstart,         &
    &                               sstart,         &
    &                               sout,           &
    &                               rout,           &
    &                               zout,           &
    &                               pcros_in
  USE cylindrical_coord_mod, ONLY : nlinp_coil_dat, &
    &                               mag_file,       &
    &                               flx_file
  USE vessel_mod,            ONLY : lvessel,        &
    &                               lcheck_vessel,  &
    &                               lvessel_vtk,    &
    &                               vessel_model,   &
    &                               vessel_file
  USE limiter_mod,           ONLY : llimiter
  USE divertor_mod,          ONLY : ldivertor
  USE plot_mod,              ONLY : fig_format,     &
    &                               fig_file,       &
    &                               contour_mode,   &
    &                               ncontb,         &
    &                               ncxy,           &
    &                               bminb,          &
    &                               bmaxb,          &
    &                               ibspec,         &
    &                               bspec

  IMPLICIT NONE

  NAMELIST /nlinp1/ lfile,          &
    &               lfigout,        &
    &               fig_format,     &
    &               fig_file,       &
    &               lpoin,          &
    &               lupdown,        &
    &               mode,           &
    &               def_start,      &
    &               plot_file_name, &
    &               lvessel,        &
    &               lcheck_vessel,  &
    &               lvessel_vtk,    &
    &               vessel_model,   &
    &               llimiter,       &
    &               ldivertor,      &
    &               lflxqnt,        &
    &               mr,             &
    &               h_in,           &
    &               lc_in,          &
    &               nstep,          &
    &               mcirc,          &
    &               rax0,           &
    &               zax0,           &
    &               drflx,          &
    &               dzflx,          &
    &               dpflx,          &
    &               rstart,         &
    &               zstart,         &
    &               pstart,         &
    &               sstart,         &
    &               sout,           &
    &               rout,           &
    &               zout,           &
    &               ntheta,         &
    &               pcros_in,       &
    &               lrmagax,        &
    &               rmagax,         &
    &               fbvadj,         &
    &               btor0,          &
    &               rbtor0,         &
    &               lcontb,         &
    &               ncontb,         &
    &               ncxy,           &
    &               bminb,          &
    &               bmaxb,          &
    &               ibspec,         &
    &               bspec
  NAMELIST /fopen/  mag_file,       &
    &               flx_file,       &
    &               vessel_file,    &
#ifdef POINXY
    &               xy_plot_file,   &
#endif
#ifdef LCOUT
    &               lc_plot_file,   &
#endif
    &               prof_file

  READ(5,nlinp1)
  !WRITE(6,nlinp1)
 
  READ(5,nlinp_coil_dat)
  !WRITE(6,nlinp_coil_dat)

  IF(lfile)THEN
    READ(5,fopen)
    !WRITE(6,fopen)
  END IF


END SUBROUTINE readin
