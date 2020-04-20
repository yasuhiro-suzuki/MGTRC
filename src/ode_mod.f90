!> @file ode_mod.f90
!------------------------------------------------------------------------------
!
! MODULE: ode_mod
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
MODULE ode_mod

  USE kind_spec

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: odeint,    &
    &       odeadm5,   &
#ifdef SSL2
    &       odedm_ssl, &
#endif
    &       rk4,       &
    &       rkf45,     &
    &       rkg6,      &
    &       dopri5

CONTAINS

  SUBROUTINE odeint (x0, y0, ln, h, lnx, method, fun, & !(in)
    &                y,                               & !(inout)
    &                icount, iout                     & !(out)
    &               )

    IMPLICIT NONE

!Arguments
    INTEGER, INTENT(IN) :: ln,  &
      &                    lnx
    INTEGER, INTENT(OUT) :: icount, &
      &                     iout
    REAL(DP), INTENT(IN) :: x0,    &
      &                     h,     &
      &                     y0(ln)
    REAL(DP), INTENT(OUT) :: y(ln,lnx)
    EXTERNAL :: method, &
      &         fun
!Local variables
    INTEGER :: i
    REAL(DP) :: x1, y1(ln)


    icount =  1
    y(:,:) =  0.0_DP
    y(:,1) =  y0(:)

    loop_main : DO i=2,lnx
      x1    =  x0 + h * (i - 2)
      y1(:) =  y(:,i-1)
      CALL method(x1, h, ln, fun, y1(:), iout)
      IF(iout == 1) EXIT
      y(:,i) =  y1(:)
      icount =  i
    END DO loop_main


    RETURN
  END SUBROUTINE odeint

  SUBROUTINE odeadm5 (x0, y0, ln, h, lnx, method, fun, & !(in)
    &                 y,                               & !(inout)
    &                 icount, iout                     & !(out)
    &                )

    IMPLICIT NONE

!Arguments
    INTEGER, INTENT(IN) :: ln,  &
      &                    lnx
    INTEGER, INTENT(OUT) :: icount, &
      &                     iout
    REAL(DP), INTENT(IN) :: x0,    &
      &                     h,     &
      &                     y0(ln)
    REAL(DP), INTENT(OUT) :: y(ln,lnx)
    EXTERNAL :: method, &
      &         fun
!Local variables
    INTEGER :: i, &
      &        j
    REAL(DP) :: error,    &
      &         x1(lnx),  &
      &         yold(ln), &
      &         y1(ln,6)


    icount =  1
    x1(1)  =  0.0_DP
    y(:,:) =  0.0_DP
    y(:,1) =  y0(:)

    loop_pre : DO i=2,5
      x1(i)   =  x0 + h * (i - 1)
      y1(:,1) =  y(:,i-1)
      CALL method(x1(i) - h, h, ln, fun, y1(:,1), iout)
      IF(iout == 1) RETURN
      y(:,i) =  y1(:,1)
      icount =  i
    END DO loop_pre

    loop_main : DO i=6,lnx
      x1(i) =  x0 + h * (i - 1)
      CALL fun(x1(i-1), y(:,i-1), y1(:,1), iout)
      CALL fun(x1(i-2), y(:,i-2), y1(:,2), iout)
      CALL fun(x1(i-3), y(:,i-3), y1(:,3), iout)
      CALL fun(x1(i-4), y(:,i-4), y1(:,4), iout)
      CALL fun(x1(i-5), y(:,i-5), y1(:,5), iout)
      !y(:,i)  =  y(:,i-1) + h * (55 * y1(:,1) - 59 * y1(:,2) + 37 * y1(:,3) - 9 * y1(:,4) ) / 24.0_DP ! 4th-order Adams–Bashforth
      y(:,i)  =  y(:,i-1) + h * (1901 * y1(:,1) - 2774 * y1(:,2) + 2616 * y1(:,3) - 1274 * y1(:,4) + 251 * y1(:,5)) / 720.0_DP ! 5th-order Adams–Bashforth
      yold(:) =  y(:,i)
      DO j=1,500
        CALL fun(x1(i), y(:,i), y1(:,6), iout)
        IF(iout == 1) RETURN
        !y(:,i) =  y(:,i-1) + h * (9 * y1(:,6) + 19 * y1(:,1) - 5 * y1(:,2) + y1(:,3)) / 24.0_DP ! 4th-order Adams–Moulton
        y(:,i) =  y(:,i-1) + h * (251 * y1(:,6) + 646 * y1(:,1) - 264 * y1(:,2) + 106 * y1(:,3) - 19 * y1(:,4)) / 720.0_DP ! 5th-order Adams–Moulton
        error  =  MAXVAL(ABS(y(:,i) - yold(:)))
        IF(error <= 1.0E-08_DP) EXIT
        yold(:) =  y(:,i)
      END DO
      icount =  i
    END DO loop_main


    RETURN
  END SUBROUTINE odeadm5

  SUBROUTINE rk4 (x0, h, n, fun, & !(in)
    &             y0,            & !(inout)
    &             iout           & !(out)
    &            )


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: x0, &
      &                     h
    REAL(DP), INTENT(INOUT) :: y0(n)
!Local variables
    REAL(DP) :: x1,      &
      &         k1(n),   &
      &         k2(n),   &
      &         k3(n),   &
      &         k4(n),   &
      &         yp(n),   &
      &         ywork(n)
 

!
!... 1-stage
!
    x1 =  x0
    CALL fun(x1, y0(:), yp(:), iout)
    IF(iout == 1) RETURN
    k1(:)    =  h * yp(:)
    ywork(:) =  y0(:) + k1(:) / 2

!
!... 2-stage
!
    x1 =  x0 + h / 2
    CALL fun(x1, ywork(:), yp(:), iout)
    IF(iout == 1) RETURN
    k2(:)    =  h * yp(:)
    ywork(:) =  y0(:) + k2(:) / 2

!
!... 3-stage
!
    x1 =  x0 + h / 2
    CALL fun(x1, ywork(:), yp(:), iout)
    IF(iout == 1) RETURN
    k3(:)    =  h * yp(:)
    ywork(:) =  y0(:) + k3(:)

!
!... 4-stage
!
    x1 =  x0 + h
    CALL fun(x1, ywork(:), yp(:), iout)
    IF(iout == 1) RETURN
    k4(:) =  h * yp(:)

    y0(:) =  y0(:) + (k1(:) + 2 * k2(:) + 2 * k3(:) + k4(:)) / 6


    RETURN 
  END SUBROUTINE rk4

  SUBROUTINE rkf45 (x0, h, n, fun, & !(in)
    &               y0,            & !(inout)
    &               iout           & !(out)
    &              )


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: x0, &
      &                     h
    REAL(DP), INTENT(INOUT) :: y0(n)
!Local variables
    INTEGER, PARAMETER :: imax =  1000
    INTEGER :: i
    REAL(DP), PARAMETER :: tol =  1.0E-12_DP
    REAL(DP) :: xend, x, x1, h1, Rmax, delta, k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), yp(n), ywork(n), R(n)


    xend =  x0 + h
    x1   =  x0
    h1   =  h

    main_loop : DO i=1,imax
!
!... 1-stage
!
      x =  x1
      CALL fun(x, y0(:), yp(:), iout)
      IF(iout == 1) RETURN
      k1(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + k1(:) / 4

!
!... 2-stage
!
      x =  x1 + h1 / 4
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k2(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 3 * k1(:) / 32 + 9 * k2(:) / 32

!
!... 3-stage
!
      x =  x1 + 3 * h1 / 8
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k3(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 1932 * k1(:) / 2197 - 7200 * k2(:) / 2197 + 7296 * k3(:) / 2197

!
!... 4-stage
!
      x =  x1 + 12 * h1 / 13
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k4(:) =  h1 * yp(:)
      ywork(:) =  y0(:) + 439 * k1(:) / 216 - 8 * k2(:) + 3680 * k3(:) / 513 - 845 * k4(:) / 4104

!
!... 5-stage
!
      x =  x1 + h1
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k5(:) =  h1 * yp(:)
      ywork(:) =  y0(:) - 8 * k1(:) / 27 + 2 * k2(:) - 3544 * k3(:) / 2565 + 1859 * k4(:) / 4104 - 11 * k5(:) / 40

!
!... 6-stage
!
      x =  x1 + h1 / 2
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k6(:) =  h1 * yp(:)

      R(:) = (k1(:) / 360 - 128 * k3(:) / 4275 - 2197 * k4(:) / 75240 + k5(:) / 50 + 2 * k6(:) / 55) / h1
      Rmax =  MAXVAL(ABS(R(:)))

      IF(Rmax < tol)THEN
        ! fifth-order
        y0(:) =  y0(:) + 16 * k1(:) / 135 + 6656 * k3(:) / 12825 + 28561 * k4(:) / 56430 - 9 * k5(:) / 50 + 2 * k6(:) / 55
        ! forth-oder
        !y0(:) =  y0(:) + 25 * k1(:) / 216 + 1408 * k3(:) / 2565 + 2197 * k4(:) / 4104 - k5(:) / 5
        x1 =  x1 + h1
        IF(x1 == xend)THEN
          EXIT main_loop
        ELSE IF(x1 > xend)THEN
          h1 =  xend - x1
        END IF
      ELSE
        delta =  0.84_DP * (tol / Rmax)**0.25_DP
        !PRINT *, i, delta
        IF(delta < 0.1_DP)THEN
          h1 =  0.1_DP * h1
        ELSE
          h1 =  delta * h1
        END IF

      END IF

    END DO main_loop


    RETURN
  END SUBROUTINE rkf45

  SUBROUTINE rkg6 (x0, h, n, fun, & !(in)
    &              y0,            & !(inout)
    &              iout           & !(out)
    &             )

!
!    This subroutine is written by Y. Suzuki
!        at Graduate School of Energy Science (Kyoto Univ)
!         2002/12/28
!
!    Based program is wrtten by K. Hamamatsu
!        at Faculty of Science (Hiroshima Univ.)
!          1980/12/18
!
!  Runge-Kutta-Huta Formulas  ( Sixth order 8-stage )
!
!  <<< Reference >>>
!  " Improved Sixth-order Runge-kutta formulas and Approximate
!    Continuous Solution of Ordinary Differential Equation "
!  by D. Sarafyan:  J. Math. Anal. Appl. 40, 436-455 (1972)
!

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: x0, &
      &                     h
    REAL(DP), INTENT(INOUT) :: y0(n)

    REAL(DP) :: x1,  &
      &         c01, &
      &         c02, &
      &         c03, &
      &         c04, &
      &         c05, &
      &         c06, &
      &         c07, &
      &         c08, &
      &         c09, &
      &         c10, &
      &         c11, &
      &         c12, &
      &         c13, &
      &         c14, &
      &         c15, &
      &         c16, &
      &         c17, &
      &         c18, &
      &         c19, &
      &         c20, &
      &         c21, &
      &         c22, &
      &         c23, &
      &         c24, &
      &         c25, &
      &         c26, &
      &         c27
    REAL(DP) :: f(n,9)
    EXTERNAL :: fun


    c01 =  h / 9.0_DP
    c02 =  h * 0.4166666666666667e-01_DP
    c03 =  h * 0.125_DP
    c04 =  h / 6.0_DP
    c05 =  h * 0.5_DP
    c06 =  h * 0.6666666666666667_DP
    c07 =  h / 3.0_DP
    c08 =  h * 0.375_DP
    c09 =  h * 0.1333333333333333e+01_DP
    c10 =  h * 0.3333333333333333e+01_DP
    c11 =  h * 0.7e+01_DP
    c12 =  h * 0.9666666666666667e+01_DP
    c13 =  h * 0.1533333333333333e+02_DP
    c14 =  h * 0.6111111111111111_DP
    c15 =  h * 0.1166666666666667e+01_DP
    c16 =  h * 0.1375e+01_DP
    c17 =  h * 0.8333333333333333_DP
    c18 =  h * 0.4390243902439024_DP
    c19 =  h * 0.8780487804878049_DP
    c20 =  h * 0.1304878048780488e+01_DP
    c21 =  h * 0.2097560975609756e+01_DP
    c22 =  h * 0.2963414634146341e+01_DP
    c23 =  h * 0.4317073170731707e+01_DP
    c24 =  h * 0.3214285714285714e-01_DP
    c25 =  h * 0.4880952380952381e-01_DP
    c26 =  h * 0.2571428571428571_DP
    c27 =  h * 0.3238095238095238_DP

    f(:,:) =  0.0_DP
    f(:,1) =  y0(:)

!
!... 1-stage
!
    x1 =  x0
    CALL fun(x1, f(:,1), f(:,3), iout)
    f(:,2) =  c01 * f(:,3) &
      &    +        f(:,1)
!
!... 2-stage
!
    x1 =  x0 + c01
    CALL fun(x1, f(:,2), f(:,4), iout)
    f(:,2) =  c02 * f(:,3) &
      &    +  c03 * f(:,4) &
      &    +        f(:,1)
!
!... 3-stage
!
    x1 =  x0 + c04
    CALL fun(x1, f(:,2), f(:,5), iout)
    f(:,2) =  c04 * f(:,3) &
      &    -  c05 * f(:,4) &
      &    +  c06 * f(:,5) &
      &    +        f(:,1)
!
!... 4-stage
!
    x1 =  x0 + c07
    CALL fun(x1, f(:,2), f(:,6), iout)
    f(:,2) =  c03 * f(:,3) &
      &    +  c08 * f(:,6) &
      &    +        f(:,1)
!
!... 5-stage
!
    x1 =  x0 + c05
    CALL fun(x1, f(:,2), f(:,7), iout)
    f(:,2) = -c09 * f(:,3) &
      &    +  c10 * f(:,7) &
      &    -  c11 * f(:,4) &
      &    -  c12 * f(:,6) &
      &    +  c13 * f(:,5) &
      &    +        f(:,1)
!
!... 6-stage
!
    x1 =  x0 + c06
    CALL fun(x1, f(:,2), f(:,8), iout)
    f(:,2) = -c01 * f(:,3) &
      &    +  c03 * f(:,8) &
      &    +  c14 * f(:,7) &
      &    -  c15 * f(:,5) &
      &    +  c16 * f(:,4) &
      &    +        f(:,1)
!
!... 7-stage
!
    x1 =  x0 + c17
    CALL fun(x1, f(:,2), f(:,9), iout)
    f(:,2) = -c18 * f(:,8) &
      &    +  c19 * f(:,9) &
      &    +  c20 * f(:,3) &
      &    -  c21 * f(:,7) &
      &    -  c22 * f(:,4) &
      &    +  c23 * f(:,6) &
      &    +        f(:,1)
!
!... 8-stage
!
    x1 =  x0 + h
    CALL fun(x1, f(:,2), f(:,4), iout)
    f(:,1) =  c24 * (f(:,6) + f(:,8)) &
      &    +  c25 * (f(:,3) + f(:,4)) &
      &    +  c26 * (f(:,5) + f(:,9)) &
      &    +  c27 *  f(:,7)           &
      &    +         f(:,1)
    y0(:)  =  f(:,1)


    RETURN
  END SUBROUTINE rkg6

  SUBROUTINE dopri5 (x0, h, n, fun, & !(in)
    &                y0,            & !(inout)
    &                iout           & !(out)
    &               )


    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: x0, &
      &                     h
    REAL(DP), INTENT(INOUT) :: y0(n)
!Local variables
    INTEGER, PARAMETER :: imax =  1000
    INTEGER :: i
    REAL(DP), PARAMETER :: tol =  1.0E-12_DP
    REAL(DP) :: xend, x, x1, h1, Rmax, delta, k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), yp(n), ywork(n), R(n)


    xend =  x0 + h
    x1   =  x0
    h1   =  h

    main_loop : DO i=1,imax
!
!... 1-stage
!
      x =  x1
      CALL fun(x, y0(:), yp(:), iout)
      IF(iout == 1) RETURN
      k1(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + k1(:) / 5

!
!... 2-stage
!
      x =  x1 + h1 / 5
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k2(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 3 * k1(:) / 40 + 9 * k2(:) / 40

!
!... 3-stage
!
      x =  x1 + 3 * h1 / 10
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k3(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 44 * k1(:) / 45 - 56 * k2(:) / 15 + 32 * k3(:) / 9

!
!... 4-stage
!
      x =  x1 + 4 * h1 / 5
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k4(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 19372 * k1(:) / 6561 - 25360 * k2(:) / 2187 + 64448 * k3(:) / 6561 - 212 * k4(:) / 729

!
!... 5-stage
!
      x =  x1 + 8 * h1 / 9
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k5(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 9017 * k1(:) / 3168 - 355 * k2(:) / 33 + 46732 * k3(:) / 5247 + 49 * k4(:) / 176 - 5103 * k5(:) / 18656

!
!... 6-stage
!
      x =  x1 + h1
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k6(:)    =  h1 * yp(:)
      ywork(:) =  y0(:) + 35 * k1(:) / 384 + 500 * k3(:) / 1113 + 125 * k4(:) / 192 - 2187 * k5(:) / 6784 + 11 * k6(:) / 84

!
!... 7-stage
!
      x =  x1 + h1
      CALL fun(x, ywork(:), yp(:), iout)
      IF(iout == 1) RETURN
      k7(:)    =  h1 * yp(:)

      R(:) = (71 * k1(:) / 57600 - 71 * k3(:) / 16695 + 71 * k4(:) / 1920 - 17253 * k5(:) / 339200 + 88 * k6(:) / 2100 - k7(:) / 40) / h1
      Rmax =  MAXVAL(ABS(R(:)))

      IF(Rmax < tol)THEN
        ! fifth-order
        y0(:) =  y0(:) + 35 * k1(:) / 384 + 500 * k3(:) / 1113 + 125 * k4(:) / 192 - 2187 * k5(:) / 6784 + 11 * k6(:) / 84
        ! forth-order
        !y0(:) =  y0(:) + 5179 * k1(:) / 57600 + 7571 * k3(:) / 16695 +393 * k4(:) / 640 - 92097 * k5(:) / 339200 + 187 * k6(:) / 2100 + k7(:) / 40

        x1 =  x1 + h1
        IF(x1 == xend)THEN
          EXIT main_loop
        ELSE IF(x1 > xend)THEN
          h1 =  xend - x1
        END IF
      ELSE
        delta =  0.84_DP * (tol / Rmax)**0.25_DP
        !PRINT *, i, Rmax, delta
        IF(delta < 0.1_DP)THEN
          h1 =  0.1_DP * h1
        ELSE
          h1 =  delta * h1
        END IF

      END IF

    END DO main_loop


    RETURN 
  END SUBROUTINE dopri5

END MODULE ode_mod
