!=trace3d.f90
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
SUBROUTINE trace2d

  USE kind_spec
  USE file_mod,              ONLY : plot_file_name
  USE param1,                ONLY : pi
  USE axis_mod,              ONLY : raxis,              &
    &                               zaxis
  USE fline_mod,             ONLY : def_start,          &
    &                               lupdown,            &
    &                               lfigout,            &
    &                               lcontb,             &
    &                               mr,                 &
    &                               nstep,              &
    &                               mcirc,              &
    &                               drflx,              &
    &                               dzflx,              &
    &                               rstart,             &
    &                               zstart,             &
    &                               pstart,             &
    &                               sstart,             &
    &                               sout,               &
    &                               rout,               &
    &                               zout,               &
    &                               pcros_in
  USE cylindrical_coord_mod, ONLY : pi2m,               &
    &                               mtor,               &
    &                               nt0b,               &
    &                               rminb,              &
    &                               rmaxb,              &
    &                               zminb,              &
    &                               zmaxb,              &
    &                               mgval1,             &
    &                               mgval3
  USE vessel_mod,            ONLY : lvessel,            &
    &                               lcheck_vessel,      &
    &                               vessel,             &
    &                               check_vessel
  USE limiter_mod,           ONLY : llimiter,           &
    &                               nphi_lim,           &
    &                               phi_lim,            &
    &                               check_limiter
  USE ode_mod,               ONLY : odeint,             &
    &                               odeadm5,            &
    &                               dopri5
  USE plot_mod,              ONLY : init_plot_poincare, &
    &                               set_plot_poincare,  &
    &                               plot_poincare,      &
    &                               init_plot_bcontour, &
    &                               end_plot_bcontour,  &
    &                               plot_bcontour,      &
    &                               init_plot_profile,  &
    &                               end_plot_profile,   &
    &                               plot_profile_reff,  &
    &                               plot_vessel
 
  IMPLICIT NONE

!Local varibales
  INTEGER, PARAMETER :: ln = 2
  INTEGER :: lnx, &
    &        i,   &
    &        j,   &
    &        k,   &
    &        l,   &
    &        js
  REAL(DP), PARAMETER :: eps =  1.0E-03_DP
  REAL(DP) :: x,        &
    &         y,        &
    &         z,        &
    &         r,        &
    &         phi,      &
    &         rq,       &
    &         zq,       &
    &         pq,       &
    &         r0,       &
    &         z0,       &
    &         p0,       &
    &         r1,       &
    &         z1,       &
    &         p1,       &
    &         p0s,      &
    &         p1s,      &
    &         w0,       &
    &         w1,       &
    &         er,       &
    &         ez,       &
    &         ep,       &
    &         ds,       &
    &         s,        &
#ifdef POINXY
    &         zcros,    &
#endif
    &         pcros(16)
  CHARACTER(LEN=100) :: fmt
  CHARACTER(LEN=300) :: plot_file
!
  REAL(DP) :: rmin, &
    &         rmax, &
    &         rmid
  REAL(DP), ALLOCATABLE :: ss(:)
!for odeint
  INTEGER, ALLOCATABLE :: icount(:,:), &
    &                     iout(:,:)
  REAL(DP) :: h
  REAL(DP), ALLOCATABLE :: x0(:),     &
    &                      f0(:,:),   &
    &                      f(:,:,:,:)
  EXTERNAL :: subf2d
!for check vessel
  INTEGER :: ivessel
!for check limiter
  INTEGER :: ilimiter


  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE trace2d                                     '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *

#ifdef POINXY
  zcros     =  0.0_DP
#endif
  pcros(:)  =  0.0_DP
  !pcros(1)  =  0.0_DP
  pcros(2)  =      pi2m / 8
  pcros(3)  =  2 * pi2m / 8
  pcros(4)  =  3 * pi2m / 8
  pcros(5)  =  4 * pi2m / 8
  pcros(6)  =  5 * pi2m / 8
  pcros(7)  =  6 * pi2m / 8
  pcros(8)  =  7 * pi2m / 8

  i =  8
  loop0010: DO l=1,8
    pcros_in(l) =  pcros_in(l) * pi / 180
    IF((mtor /= 1) .AND. (pcros_in(l) > pi2m)) pcros_in(l) =  0.0_DP
    IF(pcros_in(l) /= 0.0_DP)THEN
      i        =  i + 1
      pcros(i) =  pcros_in(l)
    END IF
  END DO loop0010

  h   =  pi2m / nstep
  lnx =  nstep * mtor * mcirc !+ 1

  ALLOCATE(x0(0:mr), f0(0:mr,ln), f(2,0:mr,ln,lnx))
  ALLOCATE(icount(2,0:mr), iout(2,0:mr))

  x0(:)       =  pi2m * pstart
  f0(:,:)     =  0.0_DP
  f(:,:,:,:)  =  0.0_DP

  icount(:,:) =  0
  iout(:,:)   =  0

  PRINT *
  PRINT *
  PRINT *, ' STEP0: CALCULATION OF STARTING POINTS'
  PRINT *
  PRINT *

  SELECT CASE(TRIM(def_start))
  
    CASE('rtheta')

      PRINT *, '  js       rho       Psi       R[m]'
      PRINT *, '--------------------------------------'

      fmt = '(I5, 2X, 3(2X,F8.5))'

      ALLOCATE(ss(0:mr))

      IF(sstart > sout)THEN
        sstart =  0.0_DP
        sout   =  1.0_DP
      END IF

      IF(sstart < 0) sstart =  0.0_DP
      IF(sout   > 1) sout   =  1.0_DP

      ds = (sout - sstart) / mr
      loop0020: DO js=0,mr
        s      =  sstart + ds * js
        ss(js) =  s**2
      END DO loop0020

      IF(ss(mr) == 1.0_DP) ss(mr) =  ss(mr) - eps

      f0(0,1) =  raxis
      f0(0,2) =  zaxis

      loop0030: DO js=1,mr

        rmin   =  raxis
        rmax   =  rmaxb

        loop0031: DO i=1,100
          rmid = (rmin + rmax) / 2

          CALL mgval3(rmid, pi2m * pstart, zaxis, s)

          IF(s < ss(js))THEN
            rmin =  rmid
          ELSE
            rmax =  rmid
          END IF
        END DO loop0031

        f0(js,1) =  rmid
        f0(js,2) =  zaxis

        !IF(js == mr) f0(js,1) =  f0(js,1) - 1.0E-04_DP

        PRINT fmt, js, SQRT(ss(js)), ss(js), f0(js,1)

      END DO loop0030

      PRINT *
      PRINT *
      PRINT *, '  mr      nstep     turn#     lnx'
      PRINT *, '-----------------------------------'

      fmt = '(I5, 3I10)'

      PRINT fmt, mr,  nstep, mcirc, lnx

      DEALLOCATE(ss)

    CASE('file', 'FILE')

      IF(ALLOCATED(x0))     DEALLOCATE(x0)
      IF(ALLOCATED(f0))     DEALLOCATE(f0)
      IF(ALLOCATED(f))      DEALLOCATE(f)
      IF(ALLOCATED(icount)) DEALLOCATE(icount)
      IF(ALLOCATED(iout))   DEALLOCATE(iout)

      READ(40,*) mr

      ALLOCATE(x0(0:mr), f0(0:mr,ln), f(2,0:mr,ln,lnx))
      ALLOCATE(icount(2,0:mr), iout(2,0:mr))

      x0(0)   =  pstart * pi2m
      f0(0,1) =  raxis
      f0(0,2) =  zaxis

      PRINT *, ' READING START POINTS FROM A FILE...'
      PRINT *
      PRINT *, '  mr      X [m]      Y [m]       Z [m]'
      PRINT *, '-----------------------------------------'

      fmt = '(I5, 3ES12.4)'

      DO js=1,mr
        READ(40,*) x, y, z
        PRINT fmt, js, x, y, z
        r        =  SQRT(x**2 + y**2)
        phi      =  ATAN2(y, x)
        x0(js)   =  phi
        f0(js,1) =  r
        f0(js,2) =  z
      END DO

      PRINT *
      PRINT *
      PRINT *, '  js      R [m]         Z [m]       phi [rad]'
      PRINT *, '-----------------------------------------------'

      fmt = '(I5, 3(2X,F12.9))'

      DO js=1,mr
        PRINT fmt, js, f0(js,1), f0(js,2), x0(js)
      END DO

      !STOP
    CASE DEFAULT

      f0(0,1) =  raxis
      f0(0,2) =  zaxis

      IF(rstart == 0.0_DP) rstart =  raxis + drflx
      IF(zstart == 0.0_DP) zstart =  zaxis + dzflx

      IF(rout /= 0.0_DP) drflx = (rout - rstart) / (mr - 1)
      IF(zout /= 0.0_DP) drflx = (zout - zstart) / (mr - 1)

      loop0040: DO js=1,mr
        f0(js,1) =  rstart + drflx * (js - 1)
        f0(js,2) =  zstart + dzflx * (js - 1)
      END DO loop0040

      PRINT *, '  Rstart [m]    Rout [m]  Zstart [m]    Zout [m]  Pstart [deg] '
      PRINT *, '---------------------------------------------------------------'

      fmt = '(5F12.4)'

      PRINT fmt, f0(1,1), f0(mr,1), f0(1,2), f0(mr,2), x0(0) * 180 / pi

      PRINT *
      PRINT *, '   dRflx [m]   dZflx [m]'
      PRINT *, '-------------------------'

      PRINT fmt, drflx, dzflx

      PRINT *
      PRINT *, '  mr      nstep     turn#     lnx'
      PRINT *, '-----------------------------------'

      fmt = '(I5, 3I10)'

      PRINT fmt, mr, nstep, mcirc, lnx

  END SELECT

  PRINT *
  PRINT *
  PRINT *, ' STEP1: FIELD LINE TRACING FROM STARTING POINTS'
  PRINT *
  PRINT *, '  mr  iout   icount   rstart   zstart '

  fmt = '(2I5, I10, 2F9.4)'

  loop0050 : DO js=0,mr

    !CALL odeint(x0(js), f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,js,:,:), icount(1,js), iout(1,js))
    CALL odeadm5(x0(js), f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,js,:,:), icount(1,js), iout(1,js))
    !IF(lupdown) CALL odeint(x0(js), f0(js,:), ln, -h, lnx, dopri5, subf2d, f(2,js,:,:), icount(2,js), iout(2,js))
    IF(lupdown) CALL odeadm5(x0(js), f0(js,:), ln, -h, lnx, dopri5, subf2d, f(2,js,:,:), icount(2,js), iout(2,js))

    PRINT fmt, js, iout(1,js), icount(1,js), f0(js,1), f0(js,2)
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

    DO js=0,mr

      r1 =  f0(js,1)
      z1 =  f0(js,2)
      p1 =  x0(js)

      IF(p1 < 0.0_DP)THEN
        k  =  ABS(p1 / pi2m) + 1
        p1 =  k * pi2m + p1
      ELSE IF(p1 > pi2m)THEN
        k  =  p1 / pi2m
        p1 =  p1 - k * pi2m
      END IF

      loop_limiter1 : DO i=2,icount(1,js)

        r0 =  r1
        z0 =  z1
        p0 =  p1

        r1 =  f(1,js,1,i)
        z1 =  f(1,js,2,i)
        p1 =  x0(js) + h * (i - 1)

        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF

        IF(lcheck_vessel)THEN
          CALL check_vessel(r1, z1, p1, ivessel)
          IF(ivessel == 0)THEN
            PRINT fmt, 1, js, i, r1, z1, p1 * 180 / pi, 'vessel'
            icount(1,js) =  i
            iout(1,js)   =  1
            EXIT loop_limiter1
          END IF
        END IF

        p0s =  p0
        p1s =  p1
        IF(ABS(p1 - p0) > 0.5_DP * pi2m)THEN
          IF(p1 < p0)THEN
            p0s =  p0 - pi2m
          ELSE
            p1s =  p1 - pi2m
          END IF
        END IF

        DO l=1,nphi_lim
          w0 =  phi_lim(l) - p0s
          w1 =  phi_lim(l) - p1s
          IF(w0 * w1 < 0.0_DP)THEN
            er = (r1 - r0) / (p1s - p0s)
            ez = (z1 - z0) / (p1s - p0s)
            ep =  phi_lim(l) - p0s
            rq =  er * ep + r0
            zq =  ez * ep + z0
            CALL check_limiter(rq, zq, l, ilimiter)
            IF(ilimiter == 1)THEN
              PRINT fmt, 1, js, i, rq, zq, phi_lim(l) * 180 / pi, 'limiter'
              icount(1,js) =  i
              iout(1,js)   =  1
              EXIT loop_limiter1
            END IF 
          END IF
        END DO
      END DO loop_limiter1
      IF(lupdown)THEN

        r1 =  f0(js,1)
        z1 =  f0(js,2)
        p1 =  x0(js)

        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF

        loop_limiter2 : DO i=2,icount(2,js)
          r0 =  r1
          z0 =  z1
          p0 =  p1

          r1 =  f(2,js,1,i)
          z1 =  f(2,js,2,i)
          p1 =  x0(js) - h * (i - 1)

          IF(p1 < 0.0_DP)THEN
            k  =  ABS(p1 / pi2m) + 1
            p1 =  k * pi2m + p1
          ELSE IF(p1 > pi2m)THEN
            k  =  p1 / pi2m
            p1 =  p1 - k * pi2m
          END IF

          IF(lcheck_vessel)THEN
            CALL check_vessel(r1, z1, p1, ivessel)
            IF(ivessel == 0)THEN
              PRINT fmt, 2, js, i, r1, z1, p1 * 180 / pi, 'vessel'
              icount(2,js) =  i
              iout(2,js)   =  1
              EXIT loop_limiter2
            END IF
          END IF

          p0s =  p0
          p1s =  p1
          IF(ABS(p1- p0) > 0.5_DP * pi2m)THEN
            IF(p1 < p0)THEN
              p0s =  p0 - pi2m
            ELSE
              p1s =  p1 - pi2m
            END IF
          END IF

          DO l=1,nphi_lim
            w0 =  phi_lim(l) - p0s
            w1 =  phi_lim(l) - p1s
            IF(w0 * w1 < 0.0_DP)THEN
              er = (r1 - r0) / (p1s - p0s)
              ez = (z1 - z0) / (p1s - p0s)
              ep =  phi_lim(l) - p0s
              rq =  er * ep + r0
              zq =  ez * ep + z0
              CALL check_limiter(rq, zq, l, ilimiter)
              IF(ilimiter == 1)THEN
                PRINT fmt, 2, js, i, rq, zq, phi_lim(l) * 180 / pi, 'limiter'
                icount(2,js) =  i
                iout(2,js)   =  1
                EXIT loop_limiter2
              END IF 
            END IF
          END DO
        END DO loop_limiter2
      END IF
    END DO

  END IF

  IF((lfigout) .AND. (lcontb))THEN
    CALL init_plot_bcontour
  END IF

  PRINT *
  PRINT *
  PRINT *, ' WRITEING PUCTURE MAPS OF FIELD LINES'
  PRINT *
  PRINT *

  PRINT *, '  l   phi[deg]'
  PRINT *, '----------------'

  fmt = '(3X, I2, 1X, F8.2)'

  loop_poloidal: DO l=1,16

    IF((l > 8) .AND. (pcros(l) == 0.0_DP)) EXIT loop_poloidal
    IF((l > 1) .AND. (nt0b == 1)) EXIT loop_poloidal

    WRITE(plot_file, fmt='(A, I0, A5)') TRIM(plot_file_name) // '.', l, '.plot'

    OPEN(50, FILE=plot_file, FORM='formatted', STATUS='unknown')

    PRINT fmt, l, pcros(l) * 180 / pi

    IF(lfigout)THEN

      IF(l <= 8)THEN
        CALL init_plot_poincare(mtor * pcros(l))
      ELSE
        CALL init_plot_poincare(pcros(l))
      END IF

      IF(lcontb)  CALL plot_bcontour(pcros(l))
      IF(lvessel) CALL plot_vessel(pcros(l))

      CALL set_plot_poincare

    END IF

    loop_radial: DO js=0,mr
      r1 =  f0(js,1)
      z1 =  f0(js,2)
      p1 =  x0(js)
      IF(p1 < 0.0_DP)THEN
        k  =  ABS(p1 / pi2m) + 1
        p1 =  k * pi2m + p1
      ELSE IF(p1 > pi2m)THEN
        k  =  p1 / pi2m
        p1 =  p1 - k * pi2m
      END IF
      !IF(p1 == pcros(l))THEN
      !  WRITE(50,'(3(ES15.7,A1))') r1, ',', z1, ',', REAL(icount(1,js)) / (mtor * nstep), ','
      !  IF(lfigout) CALL plot_poincare(r1, z1)
      !END IF
      loop0110: DO i=2,icount(1,js)
        r0 =  r1
        z0 =  z1
        p0 =  p1

        r1 =  f(1,js,1,i)
        z1 =  f(1,js,2,i)
        p1 =  x0(js) + h * (i - 1)

        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF
        IF(p1 == pcros(l))THEN
          WRITE(50,'(3(ES15.7,A1))') r1, ',', z1, ',', REAL(icount(1,js)) / (mtor * nstep), ','
          IF(lfigout) CALL plot_poincare(r1, z1)
        END IF
        p0s =  p0
        p1s =  p1
        IF(ABS(p1 - p0) > 0.5_DP * pi2m)THEN
          IF(p1 < p0)THEN
            p0s =  p0 - pi2m
          ELSE
            p1s =  p1 - pi2m
          END IF
        END IF
        w0 =  pcros(l) - p0s
        w1 =  pcros(l) - p1s
        IF(w0 * w1 < 0.0_DP)THEN
          er = (r1 - r0) / (p1s - p0s)
          ez = (z1 - z0) / (p1s - p0s)
          ep =  pcros(l) - p0s
          rq =  er * ep + r0
          zq =  ez * ep + z0
          WRITE(50,'(3(ES15.7,A1))') rq, ',', zq, ',', REAL(icount(1,js)) / (mtor * nstep), ','
          IF(lfigout) CALL plot_poincare(rq, zq)
        END IF
      END DO loop0110
      if_lupdown: IF(lupdown)THEN
        r1 =  f0(js,1)
        z1 =  f0(js,2)
        p1 =  x0(js)
        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF
        !IF(p1 == pcros(l))THEN
        !  WRITE(50,'(3(ES15.7,A1))') r1, ',', z1, ',', REAL(icount(2,js)) / (mtor * nstep), ','
        !  IF(lfigout) CALL plot_poincare(r1, z1)
        !END IF
        loop0120: DO i=2,icount(2,js)
          r0 =  r1
          z0 =  z1
          p0 =  p1

          r1 =  f(2,js,1,i)
          z1 =  f(2,js,2,i)
          p1 =  x0(js) - h * (i - 1)

          IF(p1 < 0.0_DP)THEN
            k  =  ABS(p1 / pi2m) + 1
            p1 =  k * pi2m + p1
          ELSE IF(p1 > pi2m)THEN
            k  =  p1 / pi2m
            p1 =  p1 - k * pi2m
          END IF
          IF(p1 == pcros(l))THEN
            WRITE(50,'(3(ES15.7,A1))') r1, ',', z1, ',', REAL(icount(2,js)) / (mtor * nstep), ','
            IF(lfigout) CALL plot_poincare(r1, z1)
          END IF
          p0s =  p0
          p1s =  p1
          IF(ABS(p1- p0) > 0.5_DP * pi2m)THEN
            IF(p1 < p0)THEN
              p0s =  p0 - pi2m
            ELSE
              p1s =  p1 - pi2m
            END IF
          END IF
          w0 =  pcros(l) - p0s
          w1 =  pcros(l) - p1s
          IF(w0 * w1 < 0.0_DP)THEN
            er = (r1 - r0) / (p1s - p0s)
            ez = (z1 - z0) / (p1s - p0s)
            ep =  pcros(l) - p0s
            rq =  er * ep + r0
            zq =  ez * ep + z0
            WRITE(50,'(3(ES15.7,A1))') rq, ',', zq, ',', REAL(icount(2,js)) / (mtor * nstep), ','
            IF(lfigout) CALL plot_poincare(rq, zq)
          END IF
        END DO loop0120
      END IF if_lupdown
    END DO loop_radial
    CLOSE(50)
  END DO loop_poloidal

#ifdef POINXY
  DO js=0,mr
    r1 =  f0(js,1)
    z1 =  f0(js,2)
    p1 =  x0(js)
    IF(z1 == zcros)THEN
      WRITE(51,'(2(ES15.7,A1))') r1 * COS(p1), ',', r1 * SIN(p1), ','
    END IF

    DO i=2,icount(1,js)
      r0 =  r1
      z0 =  z1
      p0 =  p1

      r1 =  f(1,js,1,i)
      z1 =  f(1,js,2,i)
      p1 =  x0(js) + h * (i - 1)

      IF(z1 == zcros)THEN
        WRITE(51,'(2(ES15.7,A1))') r1 * COS(p1), ',', r1 * SIN(p1), ','
      END IF

      w0 =  zcros - z0
      w1 =  zcros - z1
      IF(w0 * w1 < 0.0_DP)THEN
        er = (r1 - r0) / (z1 - z0)
        ep = (p1 - p0) / (z1 - z0)
        ez =  zcros - z0
        rq =  er * ez + r0
        pq =  ep * ez + p0
        WRITE(51,'(2(ES15.7,A1))') rq * COS(pq), ',', rq * SIN(pq), ','
      END IF
    END DO
    IF(lupdown)THEN
      r1 =  f0(js,1)
      z1 =  f0(js,2)
      p1 =  x0(js)

      IF(z1 == 0.0_DP)THEN
        WRITE(51,'(2(ES15.7,A1))') r1 * COS(p1), ',', r1 * SIN(p1), ','
      END IF
      DO i=2,icount(2,js)
        r0 =  r1
        z0 =  z1
        p0 =  p1

        r1 =  f(2,js,1,i)
        z1 =  f(2,js,2,i)
        p1 =  x0(js) - h * (i - 1)

        IF(z1 == zcros)THEN
          WRITE(51,'(2(ES15.7,A1))') r1 * COS(p1), ',', r1 * SIN(p1), ','
        END IF

        w0 =  zcros - z0
        w1 =  zcros - z1
        IF(w0 * w1 < 0.0_DP)THEN
          er = (r1 - r0) / (z1 - z0)
          ep = (p1 - p0) / (z1 - z0)
          ez =  zcros - z0
          rq =  er * ez + r0
          pq =  ep * ez + p0
          WRITE(51,'(2(ES15.7,A1))') rq * COS(pq), ',', rq * SIN(pq), ','
        END IF
      END DO
    END IF
  END DO
#endif

  IF((lfigout) .AND. (lcontb))THEN
    CALL end_plot_bcontour
  END IF


  DEALLOCATE(x0, f0, f)
  DEALLOCATE(icount, iout)


END SUBROUTINE trace2d
