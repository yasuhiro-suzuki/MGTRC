!> @file trace2d_rtheta.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: trace2d_rtheta
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
SUBROUTINE trace2d_rtheta

  USE kind_spec
  USE file_mod,              ONLY : plot_file_name
  USE param1,                ONLY : pi,             &
    &                               pi2
  USE axis_mod,              ONLY : raxis,          &
    &                               zaxis
  USE fline_mod,             ONLY : def_start,      &
    &                               lupdown,        &
    &                               mr,             &
    &                               nstep,          &
    &                               mcirc,          &
    &                               drflx,          &
    &                               dzflx,          &
    &                               rstart,         &
    &                               zstart,         &
    &                               pstart,         &
    &                               sstart,         &
    &                               sout,           &
    &                               rout,           &
    &                               zout
  USE cylindrical_coord_mod, ONLY : lflux,          &
    &                               pi2m,           &
    &                               mtor,           &
    &                               rmaxb,          &
    &                               mgval3
  USE vessel_mod,            ONLY : lvessel,        &
    &                               lcheck_vessel,  &
    &                               vessel,         &
    &                               check_vessel
  USE limiter_mod,           ONLY : llimiter,       &
    &                               nphi_lim,       &
    &                               phi_lim,        &
    &                               check_limiter
  USE ode_mod,               ONLY : odeint,         &
    &                               odeadm5,        &
    &                               dopri5
 
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
    &         theta
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
  REAL(DP) :: h, &
    &         x0
  REAL(DP), ALLOCATABLE :: f0(:,:),   &
    &                      f(:,:,:,:)
  EXTERNAL :: subf2d
!for check vessel
  INTEGER :: ivessel
!for check limiter
  INTEGER :: ilimiter


  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE trace2d_rtheta                              '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *

  IF(.NOT. lflux) STOP 'A parameter LFLUX must be set to TRUE in MODE:rtheta'

  x0  =  pi2m * pstart
  h   =  pi2m / nstep
  lnx =  nstep * mtor * mcirc !+ 1

  ALLOCATE(f0(0:mr,ln), f(2,0:mr,ln,lnx))
  ALLOCATE(icount(2,0:mr), iout(2,0:mr))

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

      IF(sstart <   0) sstart =  0.0_DP
      IF(sout   >   1) sout   =  1.0_DP

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

      STOP 'DEF_START:FILE mode does not work in MODE:RTHETA'

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

      PRINT fmt, f0(1,1), f0(mr,1), f0(1,2), f0(mr,2), x0 * 180 / pi

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

    !CALL odeint(x0, f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,js,:,:), icount(1,js), iout(1,js))
    CALL odeadm5(x0, f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,js,:,:), icount(1,js), iout(1,js))
    !IF(lupdown) CALL odeint(x0, f0(js,:), ln, -h, lnx, dopri5, subf2d, f(2,js,:,:), icount(2,js), iout(2,js))
    IF(lupdown) CALL odeadm5(x0, f0(js,:), ln, -h, lnx, dopri5, subf2d, f(2,js,:,:), icount(2,js), iout(2,js))

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
      p1 =  x0

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
        p1 =  x0 + h * (i - 1)

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
        p1 =  x0

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
          p1 =  x0 - h * (i - 1)

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

  PRINT *
  PRINT *
  PRINT *, ' WRITEING PUCTURE MAPS OF FIELD LINES'
  PRINT *
  PRINT *

  PRINT *, '  phi[deg]'
  PRINT *, '------------'

  fmt = '(2X, F8.2)'

  WRITE(plot_file, fmt='(A, A12)') TRIM(plot_file_name) // '.rtheta.plot'

  OPEN(50, FILE=plot_file, FORM='formatted', STATUS='unknown')

  PRINT fmt, x0 * 180 / pi

  loop_radial: DO js=0,mr
    r1 =  f0(js,1)
    z1 =  f0(js,2)
    p1 =  x0
    IF(p1 < 0.0_DP)THEN
      k  =  ABS(p1 / pi2m) + 1
      p1 =  k * pi2m + p1
    ELSE IF(p1 > pi2m)THEN
      k  =  p1 / pi2m
      p1 =  p1 - k * pi2m
    END IF
    !IF(p1 == x0)THEN
    !  CALL mgval3(r1, x0, z1, s)
    !  theta =  ATAN2(z1 - zaxis, r1 - raxis)
    !  IF(theta < 0.0_DP) theta =  theta + pi2
    !  WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(1,js)) / (mtor * nstep), ','
    !END IF
    loop0110: DO i=2,icount(1,js)
      r0 =  r1
      z0 =  z1
      p0 =  p1

      r1 =  f(1,js,1,i)
      z1 =  f(1,js,2,i)
      p1 =  x0 + h * (i - 1)

      IF(p1 < 0.0_DP)THEN
        k  =  ABS(p1 / pi2m) + 1
        p1 =  k * pi2m + p1
      ELSE IF(p1 > pi2m)THEN
        k  =  p1 / pi2m
        p1 =  p1 - k * pi2m
      END IF
      IF(p1 == x0)THEN
        CALL mgval3(r1, x0, z1, s)
        theta =  ATAN2(z1 - zaxis, r1 - raxis)
        IF(theta < 0.0_DP) theta =  theta + pi2
        WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(1,js)) / (mtor * nstep), ','
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
      w0 =  x0 - p0s
      w1 =  x0 - p1s
      IF(w0 * w1 < 0.0_DP)THEN
        er = (r1 - r0) / (p1s - p0s)
        ez = (z1 - z0) / (p1s - p0s)
        ep =  x0 - p0s
        rq =  er * ep + r0
        zq =  ez * ep + z0
        CALL mgval3(rq, x0, zq, s)
        theta =  ATAN2(zq - zaxis, rq - raxis)
        IF(theta < 0.0_DP) theta =  theta + pi2
        WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(1,js)) / (mtor * nstep), ','
      END IF
    END DO loop0110
    if_lupdown: IF(lupdown)THEN
      r1 =  f0(js,1)
      z1 =  f0(js,2)
      p1 =  x0
      IF(p1 < 0.0_DP)THEN
        k  =  ABS(p1 / pi2m) + 1
        p1 =  k * pi2m + p1
      ELSE IF(p1 > pi2m)THEN
        k  =  p1 / pi2m
        p1 =  p1 - k * pi2m
      END IF
      !IF(p1 == x0)THEN
      !  CALL mgval3(r1, x0, z1, s)
      !  theta =  ATAN2(z1 - zaxis, r1 - raxis)
      !  IF(theta < 0.0_DP) theta =  theta + pi2
      !  WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(2,js)) / (mtor * nstep), ','
      !END IF
      loop0120: DO i=2,icount(2,js)
        r0 =  r1
        z0 =  z1
        p0 =  p1

        r1 =  f(2,js,1,i)
        z1 =  f(2,js,2,i)
        p1 =  x0 - h * (i - 1)

        IF(p1 < 0.0_DP)THEN
          k  =  ABS(p1 / pi2m) + 1
          p1 =  k * pi2m + p1
        ELSE IF(p1 > pi2m)THEN
          k  =  p1 / pi2m
          p1 =  p1 - k * pi2m
        END IF
        IF(p1 == x0)THEN
          CALL mgval3(r1, x0, z1, s)
          theta =  ATAN2(z1 - zaxis, r1 - raxis)
          IF(theta < 0.0_DP) theta =  theta + pi2
          WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(2,js)) / (mtor * nstep), ','
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
        w0 =  x0 - p0s
        w1 =  x0 - p1s
        IF(w0 * w1 < 0.0_DP)THEN
          er = (r1 - r0) / (p1s - p0s)
          ez = (z1 - z0) / (p1s - p0s)
          ep =  x0 - p0s
          rq =  er * ep + r0
          zq =  ez * ep + z0
          CALL mgval3(rq, x0, zq, s)
          theta =  ATAN2(zq - zaxis, rq - raxis)
          IF(theta < 0.0_DP) theta =  theta + pi2
          WRITE(50,'(4(ES15.7,A1))') theta, ',', SQRT(s), ',', s, ',', REAL(icount(2,js)) / (mtor * nstep), ','
        END IF
      END DO loop0120
    END IF if_lupdown
  END DO loop_radial
  CLOSE(50)


  DEALLOCATE(f0, f)
  DEALLOCATE(icount, iout)


END SUBROUTINE trace2d_rtheta
