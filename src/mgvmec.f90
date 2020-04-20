!> @file mgvmec.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: mgvmec
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
SUBROUTINE mgvmec

  USE kind_spec
  USE param1,                ONLY : pi2
  USE fline_mod,             ONLY : lfigout,        &
    &                               nstep,          &
    &                               ntheta,         &
    &                               rstart,         &
    &                               zstart,         &
    &                               pstart
  USE axis_mod,              ONLY : raxis,          &
    &                               zaxis
  USE cylindrical_coord_mod, ONLY : pi2m,           &
    &                               mtor
  USE vessel_mod,            ONLY : lvessel
  USE plot_mod,              ONLY : init_plot_vmec, &
    &                               set_plot_vmec,  &
    &                               plot_vmec,      &
    &                               plot_vessel
  USE ode_mod,               ONLY : odeint, rk4

  IMPLICIT NONE

  INTEGER, PARAMETER :: ln     = 2
  INTEGER :: i,       &
    &        j,       &
    &        k,       &
    &        l,       &
    &        nstep1,  &
    &        lnx,     &
    &        icount,  &
    &        iout,    &
    &        nznt,    &
    &        nwd
  REAL(DP) :: t0,     &
    &         h,      &
    &         pcros,  &
    &         g0(ln), &
    &         g1(ln)
  REAL(DP), ALLOCATABLE :: rrs(:),      &
    &                      zrs(:),      &
    &                      phirs(:),    &
    &                      g(:,:),      &
    &                      g_axis(:,:), &
    &                      rplot(:,:),  &
    &                      zplot(:,:)
  CHARACTER(LEN=100) :: fmt
  EXTERNAL :: subf2d
!for plplot
  CHARACTER(LEN=100) :: label


  PRINT *
  PRINT *,' ----------------------------------------------------------------'
  PRINT *,'          SUBROUTINE mgvmec                                      '
  PRINT *,' ----------------------------------------------------------------'
  PRINT *

  nstep1 =  nstep + 1
  nwd    =  nstep * mtor + 1
  lnx    =  nstep * ntheta

  h      =  pi2m / nstep

  ALLOCATE(g(ln,lnx), g_axis(ln,nwd), rrs(nwd), zrs(nwd), phirs(nwd), rplot(ntheta,nstep1), zplot(ntheta,nstep1))
 
  t0     =  pi2m * pstart
  g0(1)  =  raxis
  g0(2)  =  zaxis

  CALL odeint(t0, g0, ln, h, nwd, rk4, subf2d, g_axis, icount, iout)
  IF(nwd /= icount) STOP ' nwd /= icount'

  loop010 : do i=1,nwd
    rrs(i)   =  g_axis(1,i)
    zrs(i)   =  g_axis(2,i)
    phirs(i) =  h * (i - 1)
  END do loop010

  IF(rstart == 0.0_DP) rstart =  raxis
  IF(zstart == 0.0_DP) zstart =  zaxis

  t0     =  pi2m * pstart
  g0(1)  =  rstart
  g0(2)  =  zstart

  CALL odeint(t0, g0, ln, h, lnx, rk4, subf2d, g, icount, iout)
  IF(lnx /= icount) STOP ' lnx /= icount'

  loop110 : do l=1,nstep

    pcros =  h * (l - 1)

    loop120 : do i=1,ntheta
      k          = (i - 1) * nstep + l
      rplot(i,l) =  g(1,k)
      zplot(i,l) =  g(2,k)
    END do loop120

    IF(lfigout)THEN
      CALL init_plot_vmec(pcros)
      CALL set_plot_vmec
      CALL plot_vmec(ntheta, rplot(:,l), zplot(:,l), rrs(l), zrs(l))
      IF(lvessel) CALL plot_vessel(pcros)
    END IF

  END do loop110

  
  rplot(:,nstep1) =  rplot(:,1)
  zplot(:,nstep1) =  zplot(:,1)

  fmt = '(A8, 2F7.4)'

  PRINT *
  PRINT *
  PRINT fmt,' hout = ', MAXVAL(rplot(:,1))
  PRINT fmt,' hin  = ', MINVAL(rplot(:,1))
  PRINT fmt,' vout = ', MAXVAL(rplot(:,nstep/2+1))
  PRINT fmt,' vin  = ', MINVAL(rplot(:,nstep/2+1))
  PRINT *
  PRINT *

  nznt =  nstep * ntheta

  WRITE(88) mtor, ntheta, nstep, nznt, nwd
  WRITE(88) ((rplot(i,j), i=1,ntheta), j=1,nstep), &
    &       ((zplot(i,j), i=1,ntheta), j=1,nstep)
  WRITE(88) (phirs(i), i=1,nwd), &
    &       (rrs(i),   i=1,nwd), &
    &       (zrs(i),   i=1,nwd)

  DEALLOCATE(g, g_axis, rrs, zrs, phirs, rplot, zplot)


END SUBROUTINE mgvmec
