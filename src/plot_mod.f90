!=plot_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! SUBROUTINES for plotting
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
! These routines are coded based on the PLPLOT Graphic Libarry
! see plplot.sourceforge.net
!
!==TODO
!
! None
!

MODULE plot_mod

  USE kind_spec
  USE param1,                ONLY : pi
  USE cylindrical_coord_mod, ONLY : rminb, &
    &                               rmaxb, &
    &                               zminb, &
    &                               zmaxb, &
    &                               mgval1
  USE vessel_mod,            ONLY : nu,    &
    &                               vessel
#ifdef PLPLOT
  USE plplot
#endif

#ifdef PLPLOT77
#define PLPLOT
#endif

  IMPLICIT NONE

  PRIVATE

!
!
  INTEGER :: ncontb  =  0,  &
    &        ncxy    =  51, &
    &        ibspec  =  0
  REAL(DP) :: bmaxb       =  5.0_DP, &
    &         bminb       =  1.0_DP, &
    &         bspec(10)   =  0.0_DP
  CHARACTER(LEN=20) :: contour_mode =  'surface'
  CHARACTER(LEN=300) :: fig_format  =  'ps',        &
    &                   fig_file    =  'default.ps'
#ifdef PLPLOT
!for plplot
#if defined PLPLOT77
  INTEGER :: mag_color,  &
    &        base_cont,  &
    &        fill_width, &
    &        cont_color, &
    &        cont_width, &
    &        opt(2)
#else
  INTEGER :: cont_color
  REAL(DP) :: fill_width, &
    &         cont_width
#endif
  REAL(DP) :: xmin,  &
    &         xmax,  &
    &         ymin,  &
    &         ymax,  &
    &         br,    &
    &         bz,    &
    &         bt,    &
    &         drb,   &
    &         dzb,   &
    &         dcont, &
    &         tr(6)
  REAL(DP), ALLOCATABLE :: xx(:),     &
    &                      yy(:),     &
    &                      clevel(:), &
    &                      rgb(:),    &
    &                      zgb(:),    &
    &                      bcont(:,:)
  CHARACTER(LEN=1) :: defined
  CHARACTER(LEN=100) :: label
  COMMON /plplot/ tr
#endif

  PUBLIC :: fig_format,         &
    &       fig_file,           & 
    &       contour_mode,       &
    &       ncontb,             &
    &       ncxy,               &
    &       bminb,              &
    &       bmaxb,              &
    &       ibspec,             &
    &       bspec,              &
    &       init_plot,          &
    &       end_plot,           &
    &       init_plot_poincare, &
    &       set_plot_poincare,  &
    &       plot_poincare,      &
    &       init_plot_bcontour, &
    &       end_plot_bcontour,  &
    &       plot_bcontour,      &
    &       init_plot_profile,  &
    &       end_plot_profile,   &
    &       plot_profile_reff,  &
    &       plot_profile_R,     &
    &       plot_profile_phi,   &
    &       plot_profile_psi,   &
    &       init_plot_vmec,     &
    &       set_plot_vmec,      &
    &       plot_vmec,          &
    &       plot_vessel

CONTAINS

  SUBROUTINE init_plot

    IMPLICIT NONE


#ifdef PLPLOT
    CALL plsdev(fig_format)
    CALL plsfnam(fig_file)

    !CALL plspage(600.0_DP, 600.0_DP, 1189, 841, 10, 10)

    CALL plssub(2, 2)

    !CALL plscol0(0, 190, 190, 190)
    !CALL plscol0(7, 0, 0, 0)

    CALL plscol0(0, 255, 255, 255)
    CALL plscol0(1, 0, 0, 0)
    CALL plscol0(15, 255, 0, 0)

    CALL plinit()

    CALL plfont(2)
#endif


  END SUBROUTINE init_plot

  SUBROUTINE end_plot

    IMPLICIT NONE


#ifdef PLPLOT
    CALL plend()
#endif


  END SUBROUTINE end_plot

  SUBROUTINE init_plot_poincare (pcros & !(in)
    &                           )


    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: pcros
!Local variables
#ifdef PLPLOT
    CHARACTER(LEN=100) :: label
#endif


#ifdef PLPLOT
    CALL plcol0(1)
    CALL pLENv(rminb, rmaxb, zminb, zmaxb, 1, 2)

    WRITE(label, fmt = '(A8, F7.3)') ' M#gf = ', 180.0_DP *  pcros / pi

    CALL pllab('R [m]', 'Z [m]', label)
#endif
 

  END SUBROUTINE init_plot_poincare
  
  SUBROUTINE set_plot_poincare

    IMPLICIT NONE


#ifdef PLPLOT
    CALL plcol0(15)

    CALL plssym(0.0_DP, 0.3_DP)
#endif


  END SUBROUTINE set_plot_poincare

  SUBROUTINE plot_poincare (r, z & !(in)
    &                      )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: r, &
      &                     z
#if !defined PLPLOT77
    REAL(DP) :: r1(1), &  !local to satisfy new plplot interface requirements
      &         z1(1)     !and to keep interface for plot_poincare-call stable.
#endif


#ifdef PLPLOT
#if defined PLPLOT77
    CALL plpoin(1, r, z, 17)
#else
    r1(1) =  r
    z1(1) =  z
    CALL plpoin(r1, z1, 17)
#endif
#endif

 
  END SUBROUTINE plot_poincare

  SUBROUTINE init_plot_bcontour

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
 

#ifdef PLPLOT
    IF(ncontb <= 0) ncontb =  31
    IF(ncxy   <= 0) ncxy   =  51

    ALLOCATE(clevel(ncontb), rgb(ncxy), zgb(ncxy), bcont(ncxy,ncxy))

    drb = (rmaxb - rminb) / (ncxy - 1)
    dzb = (zmaxb - zminb) / (ncxy - 1)

    DO i=1,ncxy
      rgb(i) =  rminb + drb * (i - 1)
    END DO

    DO j=1,ncxy
      zgb(j) =  zminb + dzb * (j - 1)
    END DO

    tr(1) =  drb
    tr(2) =  0.0_DP
    tr(3) =  rminb
    tr(4) =  0.0_DP
    tr(5) =  dzb
    tr(6) =  zminb
#endif
 
    
  END SUBROUTINE init_plot_bcontour

  SUBROUTINE plot_bcontour (pcros & !(in)
    &                      )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: pcros
!Local variables
    INTEGER :: i, &
      &        j


#ifdef PLPLOT
    IF(bminb == 0.0_DP) bminb =  0.0_DP
    IF(bmaxb == 0.0_DP) bmaxb =  10.0_DP

    DO j=1,ncxy
      DO i=1,ncxy
        CALL mgval1(rgb(i), pcros, zgb(j), br, bt, bz, bcont(i,j))
        IF(bcont(i,j) < bminb) bcont(i,j) =  bminb
        IF(bcont(i,j) > bmaxb) bcont(i,j) =  bmaxb
      END DO
    END DO

    dcont = (bmaxb - bminb) / (ncontb - 1)

    DO i=1,ncontb
      clevel(i) = bminb + dcont * (i - 1)
    END DO

    SELECT CASE(TRIM(contour_mode))

      CASE('lines')

        CALL plcol0(5)
#if defined PLPLOT77
        CALL plcont(bcont(:,:), ncxy, ncxy, 1, ncxy, 1, ncxy, clevel(:), ncontb)
#else
        CALL plcont(bcont(:,:), 1, ncxy, 1, ncxy, clevel(:))
#endif

      CASE('surface')

        fill_width =  2
        cont_color =  0
        cont_width =  0

#if defined PLPLOT77
        CALL plshades0(bcont(:,:), ncxy, ncxy, defined, rminb, rmaxb, zminb, zmaxb, clevel(:), ncontb, fill_width, cont_color, cont_width, ncxy)
#else
        CALL plshades(bcont(:,:), rminb, rmaxb, zminb, zmaxb, clevel(:), fill_width, cont_color, cont_width, .true.)
#endif

      CASE DEFAULT

        CALL plcol0(5)
#if defined PLPLOT77
        CALL plcont(bcont(:,:), ncxy, ncxy, 1, ncxy, 1, ncxy, clevel(:), ncontb)
#else
        CALL plcont(bcont(:,:), 1, ncxy, 1, ncxy, clevel(:))
#endif

        fill_width =  2
        cont_color =  0
        cont_width =  0

#if defined PLPLOT77
        CALL plshades0(bcont(:,:), ncxy, ncxy, defined, rminb, rmaxb, zminb, zmaxb, clevel(:), ncontb, fill_width, cont_color, cont_width, ncxy)
#else
        CALL plshades(bcont(:,:), rminb, rmaxb, zminb, zmaxb, clevel(:), fill_width, cont_color, cont_width, .true.)
#endif

    END SELECT

    IF(ibspec > 0)THEN
      loop6040 : DO i=1,ibspec
        CALL plcol0(4)
#if defined PLPLOT77
        CALL plcont(bcont(:,:), ncxy, ncxy, 1, ncxy, 1, ncxy, bspec(i), 1)
#else
        CALL plcont(bcont(:,:), 1, ncxy, 1, ncxy, bspec(i:i))
#endif
      END DO loop6040
    END IF
#endif


  END SUBROUTINE plot_bcontour

  SUBROUTINE end_plot_bcontour

    IMPLICIT NONE


#ifdef PLPLOT
    DEALLOCATE(clevel, rgb, zgb, bcont)
#endif


  END SUBROUTINE end_plot_bcontour

  SUBROUTINE init_plot_profile (no & !(in)
    &                          )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: no


#ifdef PLPLOT
    CALL plssym(0.0_DP, 1.0_DP)

    ALLOCATE(xx(0:no), yy(0:no))
#endif


  END SUBROUTINE init_plot_profile

  SUBROUTINE plot_profile_reff (no, rav1, rav2, iota1, vp, volume, lc & !(in)
    &                          )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: no
    REAL(DP), INTENT(IN) :: rav1(0:no),   &
      &                     rav2(0:no),   &
      &                     iota1(0:no),  &
      &                     vp(0:no),     &
      &                     volume(0:no), &
      &                     lc(0:no)
!Local variable
    INTEGER :: i


#ifdef PLPLOT
    xx(:) =  rav1(0:no)

    xmin  =  0.0_DP
    xmax  =  MAXVAL(rav2(0:no)) !* 1.2_DP

! plot of iota vs r_eff

    yy(:) =  iota1(0:no)

    ymin  =  MINVAL(yy(:)) - 0.1_DP
    ymax  =  MAXVAL(yy(:)) + 0.1_DP

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '#gi / 2#gp', 'rotational transform')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of q vs r_eff

    DO i=0,no
      IF(ABS(iota1(i)) /= 0.0_DP)THEN
        yy(i) =  1.0_DP / ABS(iota1(i))
      ELSE
        yy(i) =  10.0_DP
      END IF
    END DO

    ymin  =  0.0_DP
    ymax  =  MAXVAL(yy(:)) + 1

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'q', 'safty factor')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of well vs r_eff

    yy(:) =  (vp(0) - vp(0:no)) / vp(0)

    ymin  =  MINVAL(yy(:))
    ymax  =  MAXVAL(yy(:))

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '(V#dp#u(0)-V#dp#u(a))/V#dp#u(0)', 'Magnetic Well Depth')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of volume vs r_eff

    yy(:) =  volume(0:no)

    ymin  =  0.0_DP
    ymax  =  MAXVAL(yy(:)) + 1

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'V', 'Plasma Volume')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of lc vs r_eff

    yy(:) =  lc(0:no)

    ymin  =  MINVAL(yy(:))
    ymax  =  MAXVAL(yy(:))

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'L#dC', 'Connection Length')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif
#endif


  END SUBROUTINE plot_profile_reff

  SUBROUTINE plot_profile_R (no, R, iota1, vp & !(in)
    &                       )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: no
    REAL(DP), INTENT(IN) :: R(0:no),     &
      &                     iota1(0:no), &
      &                     vp(0:no)


#ifdef PLPLOT
    xx(:) =  R(0:no)

    xmin  =  0.0_DP
    xmax  =  MAXVAL(R(0:no)) !* 1.2_DP

! plot of iota vs R

    yy(:) =  iota1(0:no)

    ymin  =  MINVAL(yy(:)) - 0.1_DP
    ymax  =  MAXVAL(yy(:)) + 0.1_DP

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '#gi / 2#gp', 'rotational transform')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of q vs R

    yy(:) =  1.0_DP / ABS(iota1(0:no))

    ymin  =  0.0_DP
    ymax  =  MAXVAL(yy(:)) + 1

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'q', 'safty factor')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of well vs R

    yy(:) =  (vp(0) - vp(0:no)) / vp(0)

    ymin  =  MINVAL(yy(:))
    ymax  =  MAXVAL(yy(:))

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '(V#dp#u(0)-V#dp#u(a))/V#dp#u(0)', 'Magnetic Well Depth')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif
#endif


  END SUBROUTINE plot_profile_R
 
  SUBROUTINE plot_profile_phi (no, phit, iota1, vp & !(in)
    &                         )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: no
    REAL(DP), INTENT(IN) :: phit(0:no),  &
      &                     iota1(0:no), &
      &                     vp(0:no)


#ifdef PLPLOT
    xx(:) =  phit(0:no) / phit(no)

    xmin  =  0.0_DP
    xmax  =  1.0_DP

! plot of iota vs s

    yy(:) =  iota1(0:no)

    ymin  =  MINVAL(yy(:)) - 0.1_DP
    ymax  =  MAXVAL(yy(:)) + 0.1_DP

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '#gi / 2#gp', 'rotational transform')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of q vs s

    yy(:) =  1.0_DP / ABS(iota1(0:no))

    ymin  =  0.0_DP
    ymax  =  MAXVAL(yy(:)) + 1

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'q', 'safty factor')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of well vs s

    yy(:) =  (vp(0) - vp(0:no)) / vp(0)

    ymin  =  MINVAL(yy(:))
    ymax  =  MAXVAL(yy(:))

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '(V#dp#u(0)-V#dp#u(a))/V#dp#u(0)', 'Magnetic Well Depth')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif
#endif


  END SUBROUTINE plot_profile_phi

  SUBROUTINE plot_profile_psi (no, psip, iota1, vp & !(in)
    &                         )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: no
    REAL(DP), INTENT(IN) :: psip(0:no),  &
      &                     iota1(0:no), &
      &                     vp(0:no)


#ifdef PLPLOT
    xx(:) =  psip(0:no) / psip(no)

    xmin  =  0.0_DP
    xmax  =  1.0_DP

! plot of iota vs Psi

    yy(:) =  iota1(0:no)

    ymin  =  MINVAL(yy(:)) - 0.1_DP
    ymax  =  MAXVAL(yy(:)) + 0.1_DP

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '#gi / 2#gp', 'rotational transform')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of q vs Psi

    yy(:) =  1.0_DP / ABS(iota1(0:no))

    ymin  =  0.0_DP
    ymax  =  MAXVAL(yy(:)) + 1

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', 'q', 'safty factor')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif

! plot of well vs Psi

    yy(:) =  (vp(0) - vp(0:no)) / vp(0)

    ymin  =  MINVAL(yy(:))
    ymax  =  MAXVAL(yy(:))

#if defined PLPLOT77
    CALL plwid(1)
#else
    CALL plwidth(1.0)
#endif
    CALL plcol0(1)
    CALL plenv(xmin, xmax, ymin, ymax, 0, 2)

    CALL pllab('r#deff#u (m)', '(V#dp#u(0)-V#dp#u(a))/V#dp#u(0)', 'Magnetic Well Depth')

    CALL plcol0(15)
#if defined PLPLOT77
    CALL plwid(3)
    CALL plline(no+1, xx(0:no), yy(0:no))
    CALL plpoin(no+1, xx(0:no), yy(0:no), 9)
#else
    CALL plwidth(3.0)
    CALL plline(xx(0:no), yy(0:no))
    CALL plpoin(xx(0:no), yy(0:no), 9)
#endif
#endif


  END SUBROUTINE plot_profile_psi

  SUBROUTINE end_plot_profile

    IMPLICIT NONE


#ifdef PLPLOT
    DEALLOCATE(xx, yy)
#endif


  END SUBROUTINE end_plot_profile

  SUBROUTINE init_plot_vmec(pcros & !(in)
    &                      )


    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: pcros
!Local variables
    CHARACTER(LEN=100) :: label


#ifdef PLPLOT
    CALL plcol0(1)
    CALL pLENv(rminb, rmaxb, zminb, zmaxb, 1, 2)

    WRITE(label, fmt = '(A8, F7.3)') ' M#gf = ', 180.0_DP *  pcros / pi

    CALL pllab('R [m]', 'Z [m]', label)
#endif


  END SUBROUTINE init_plot_vmec

  SUBROUTINE set_plot_vmec

    IMPLICIT NONE


#ifdef PLPLOT
    CALL plcol0(15)

    CALL plssym(0.0_DP, 0.3_DP)
#endif


  END SUBROUTINE set_plot_vmec

  SUBROUTINE plot_vmec(nplot, rplot, zplot, r0, z0 & !(in)
    &                 )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: nplot
    REAL(DP), INTENT(IN) :: r0,           &
      &                     z0,           &
      &                     rplot(nplot), & 
      &                     zplot(nplot)
#if !defined PLPLOT77
    REAL(DP) :: r1(1), & ! local to satisfy new plplot interface requirements
      &         z1(1)    ! and to keep interface for plot_poincare-call stable.
#endif


#ifdef PLPLOT
    CALL plcol0(2)
    CALL plcol0(15)

#if defined PLPLOT77
    CALL plpoin(nplot, rplot(:), zplot(:), 1)
#else
    CALL plpoin(rplot(:), zplot(:), 1)
#endif

    CALL plcol0(3)
#if defined PLPLOT77
    CALL plpoin(1, r0, z0, 2)
#else
    r1(1) =  r0
    z1(1) =  z0
    CALL plpoin(r1, z1, 2)
#endif
#endif


  END SUBROUTINE plot_vmec

  SUBROUTINE plot_vessel (pcros & !(in)
    &                    )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: pcros
!Local variables
    REAL(DP) :: rv(nu), &
      &         zv(nu)


#ifdef PLPLOT
    CALL vessel(pcros, rv(:), zv(:))

    CALL plcol0(3)

#if defined PLPLOT77
    CALL plline(nu, rv(:), zv(:))
#else
    CALL plline(rv(:), zv(:))
#endif
#endif


  END SUBROUTINE plot_vessel

END MODULE plot_mod
