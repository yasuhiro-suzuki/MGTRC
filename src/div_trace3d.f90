!> @file div_trace3d.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: div_trace3d
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
SUBROUTINE div_trace3d

  USE kind_spec
  USE param1,                ONLY : pi
  USE fline_mod,             ONLY : mr,            &
    &                               h_in,          &
    &                               lc_in
  USE cylindrical_coord_mod, ONLY : lflux,         &
    &                               pi2m,          &
    &                               mgval3
  USE vessel_mod,            ONLY : lvessel,       &
    &                               lcheck_vessel, &
    &                               vessel,        &
    &                               check_vessel
  USE limiter_mod,           ONLY : llimiter,      &
    &                               nphi_lim,      &
    &                               phi_lim,       &
    &                               check_limiter
  USE ode_mod,               ONLY : odeint,        &
    &                               odeadm5,       &
    &                               dopri5
 
  IMPLICIT NONE

!Local varibales
  INTEGER, PARAMETER :: ln = 3
  INTEGER :: lnx, &
    &        i,   &
    &        k,   &
    &        js
  REAL(DP) :: r,           &
    &         phi,         &
    &         z,           &
    &         r0,          &
    &         z0,          &
    &         p0,          &
    &         r1,          &
    &         z1,          &
    &         p1,          &
    &         s,           &
    &         smin
  REAL(DP), ALLOCATABLE :: xx(:), &
    &                      yy(:), &
    &                      zz(:)
  CHARACTER(LEN=100) :: fmt
!
!for odeint
  INTEGER, ALLOCATABLE :: icount(:), &
    &                     iout(:)
  REAL(DP) :: h, &
    &         x0
  REAL(DP), ALLOCATABLE :: f0(:,:),   &
    &                      f(:,:,:)
  EXTERNAL :: subf3d
!for check vessel
  INTEGER :: ivessel
!for check limiter
  INTEGER :: ilimiter


  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE div_trace3d                                 '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *


  READ(40,*) mr

  x0  =  0.0_DP
  h   =  h_in 
  lnx =  lc_in / h_in + 1

  ALLOCATE(f0(mr,ln), f(mr,ln,lnx))
  ALLOCATE(icount(mr), iout(mr))
  ALLOCATE(xx(mr), yy(mr), zz(mr))

  f0(:,:)     =  0.0_DP
  f(:,:,:)  =  0.0_DP

  icount(:) =  0
  iout(:)   =  0

  PRINT *
  PRINT *
  PRINT *, ' STEP0: CALCULATION OF STARTING POINTS'
  PRINT *
  PRINT *

  PRINT *, ' READING START POINTS FROM A FILE...'
  PRINT *
  PRINT *, '  mr      X [m]      Y [m]       Z [m]'
  PRINT *, '-----------------------------------------'

  fmt = '(I5, 3ES12.4)'

  loop0040: DO js=1,mr
    READ(40,*) xx(js), yy(js), zz(js)
    PRINT fmt, js, xx(js), yy(js), zz(js)
    r        =  SQRT(xx(js)**2 + yy(js)**2)
    phi      =  ATAN2(yy(js), xx(js))
    f0(js,1) =  r
    f0(js,2) =  zz(js)
    f0(js,3) =  phi
  END DO loop0040

  PRINT *
  PRINT *
  PRINT *, '  js      R [m]         Z [m]       phi [rad]'
  PRINT *, '-----------------------------------------------'

  fmt = '(I5, 3(2X,F12.9))'

  DO js=1,mr
    PRINT fmt, js, f0(js,1), f0(js,2), f0(js,3)
  END DO
 
  PRINT *
  PRINT *, '  mr      h [m]      Lc [m]      lnx'
  PRINT *, '----------------------------------------'

  fmt = '(I5, 2ES12.4, I10)'

  PRINT fmt, mr,  h_in, lc_in, lnx
 
  PRINT *
  PRINT *
  PRINT *, ' STEP1: FIELD LINE TRACING FROM STARTING POINTS'
  PRINT *
  PRINT *, '  mr  iout   icount   rstart   zstart '

  fmt = '(2I5, I10, 2F9.4)'

  loop0050 : DO js=1,mr

    !CALL odeint(x0, f0(js,:), ln, -h, lnx, dopri5, subf3d, f(js,:,:), icount(js), iout(js))
    CALL odeadm5(x0, f0(js,:), ln, -h, lnx, dopri5, subf3d, f(js,:,:), icount(js), iout(js))

    PRINT fmt, js, iout(js), icount(js), f0(js,1), f0(js,2)
  END DO loop0050

  IF(llimiter .OR. lcheck_vessel)THEN

    PRINT *
    PRINT *
    PRINT *, ' STEP2: CHECK HITTING POINT ON LIMITER OR VESSEL'
    PRINT *
    PRINT *

    PRINT *, '  parity   js    icount      R[m]     Z[m]    phi[deg]'
    PRINT *, '------------------------------------------------------------------'

    fmt = '(3X, I5, 1X, I5, I10, 2X, 3(2X,F7.3), 5X, A)'

    DO js=1,mr

      s    = -10.0_DP
      smin =  10.0_DP

      r1   =  f0(js,1)
      z1   =  f0(js,2)
      p1   =  f0(js,3)

      IF(p1 < 0.0_DP)THEN
        k  =  ABS(p1 / pi2m) + 1
        p1 =  k * pi2m + p1
      ELSE IF(p1 > pi2m)THEN
        k  =  p1 / pi2m
        p1 =  p1 - k * pi2m
      END IF

      IF(lflux) CALL mgval3(r1, p1, z1, s)
      IF(smin > s) smin =  s

      loop_limiter1 : DO i=2,icount(js)

        r0 =  r1
        z0 =  z1
        p0 =  p1

        r1 =  f(js,1,i)
        z1 =  f(js,2,i)
        p1 =  f(js,3,i)

        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF

        IF(lflux) CALL mgval3(r1, p1, z1, s)
        IF(smin > s) smin =  s

        IF(lcheck_vessel)THEN
          CALL check_vessel(r1, z1, p1, ivessel)
          IF(ivessel == 0)THEN
            PRINT fmt, 1, js, i, r1, z1, p1 * 180 / pi, 'vessel'
            icount(js) =  i
            iout(js)   =  1
            EXIT loop_limiter1
          END IF
        END IF

      END DO loop_limiter1

      WRITE(41,'(I5,A,3(F16.10,A))') js, ',', xx(js), ',', h * icount(js), ',', smin, ','

    END DO

  END IF


  DEALLOCATE(f0, f)
  DEALLOCATE(icount, iout)
  DEALLOCATE(xx, yy, zz)


END SUBROUTINE div_trace3d
