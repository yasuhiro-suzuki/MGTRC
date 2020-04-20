!> @file mgaxis.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: mgaxis
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
SUBROUTINE mgaxis

  USE kind_spec
  USE param1,                ONLY : pi2
  USE axis_mod,              ONLY : rax0,   &
    &                               zax0,   &
    &                               raxis,  &
    &                               zaxis,  &
    &                               baxis,  &
    &                               bpaxis, &
    &                               rax_av, &
    &                               zax_av, &
    &                               bax_av, &
    &                               bpax_av
  USE fline_mod,             ONLY : pstart
  USE cylindrical_coord_mod, ONLY : mtor,   &
    &                               pi2m,   & 
    &                               mgval1
  USE ode_mod,               ONLY : odeint, &
    &                               dopri5

  IMPLICIT NONE

!Local varibales
  INTEGER, PARAMETER :: ln    =  2,   &
    &                   lstep =  360, &
    &                   mcirc =  100
  INTEGER :: i,       &
    &        j,       &
    &        k,       &
    &        lnx,     &
    &        icount,  &
    &        iout,    &
    &        nplt,    &
    &        icros(8)
  REAL(DP) :: raxo,    &
    &         zaxo,    &
    &         rmid,    &
    &         zmid,    &
    &         eps,     &
    &         h,       &
    &         x0,      &
    &         r,       &
    &         z,       &
    &         phi,     &
    &         dr,      &
    &         dz,      &
    &         dphi,    &
    &         br,      &
    &         bp,      &
    &         bz,      &
    &         bb,      &
    &         s,       &
    &         pp
  REAL(DP), ALLOCATABLE :: f0(:), &
    &                      f(:,:)
  CHARACTER(LEN=100) :: fmt
  EXTERNAL :: subf2d


  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE mgaxis                                      '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *

  lnx   =  lstep * mcirc

  ALLOCATE(f0(ln), f(ln,lnx))

  h     =  pi2  / lstep
  x0    =  pi2m * pstart

  fmt =  '(A10, I12)'

  PRINT fmt, ' lstep = ', lstep
  PRINT fmt, ' mcirc = ', mcirc
  PRINT fmt, ' lnx   = ', lnx

  fmt =  '(A10, F12.4)'
  PRINT fmt, ' h     = ', h
  PRINT *
  PRINT *

  f0(1)  =  rax0
  f0(2)  =  zax0

  raxis  =  rax0
  zaxis  =  zax0

  PRINT *, '  searching magnetic axis'
  PRINT *
  PRINT *, ' m       r           z       phi [deg]'

  fmt =  '(I3, 2ES12.4, F10.2)'

  loop_010 : DO i=1,100
    PRINT fmt, i, raxis, zaxis, 360.0_DP * x0 / pi2
    raxo  =  raxis
    zaxo  =  zaxis
    rmid  =  raxis
    zmid  =  zaxis
    f0(1) =  raxis
    f0(2) =  zaxis
    CALL odeint(x0, f0, ln, h, lnx, dopri5, subf2d, f, icount, iout)
    IF(lnx /= icount) STOP ' lnx /= icount'
    loop020 : DO j=1,lnx
      nplt =  MOD(j, lstep)
      IF(nplt == 1)THEN
        rmid =  rmid + f(1,j)
        zmid =  zmid + f(2,j)
      END IF
    END DO loop020
    raxis =  rmid / (mcirc + 1)
    zaxis =  zmid / (mcirc + 1)
    eps   =  SQRT((raxis - raxo)**2 + (zaxis - zaxo)**2)
    IF(eps < 1.0E-06_DP)THEN
      PRINT *
      PRINT *, 'magnetic axis searching is convergenced'
      PRINT *
      EXIT loop_010
    END IF
  END DO loop_010

  CALL mgval1(raxis, x0, zaxis, br, bpaxis, bz, baxis)

  fmt =  '(//3(A7, F9.4), A8, F9.4//)'
  PRINT fmt, '  R  = ', raxis, '  Z  = ', zaxis, '  B  = ', baxis, '  Bp  = ', bpaxis

  icros(1) =  1
  icros(2) =      (lstep / mtor) / 8 + 1
  icros(3) =  2 * (lstep / mtor) / 8 + 1
  icros(4) =  3 * (lstep / mtor) / 8 + 1
  icros(5) =  4 * (lstep / mtor) / 8 + 1
  icros(6) =  5 * (lstep / mtor) / 8 + 1
  icros(7) =  6 * (lstep / mtor) / 8 + 1
  icros(8) =  7 * (lstep / mtor) / 8 + 1

  lnx   =  lstep + 1

  f0(1) =  raxis
  f0(2) =  zaxis

  CALL odeint(x0, f0, ln, h, lnx, dopri5, subf2d, f, icount, iout)

  PRINT *, '   k  phi[deg]    R[m]     Z[m]     B[T]     Bp[T]'
  PRINT *, '----------------------------------------------------'

  fmt =  '(I5, F9.2, 1X, 4F9.4)'
  loop110 : DO k=1,8
    i   =  icros(k)
    phi =  x0 + h * (i - 1)
    CALL mgval1(f(1,i), phi, f(2,i), br, bp, bz, bb)
    PRINT fmt, i, phi * 360 / pi2, f(1,i), f(2,i), bb, bp
  END DO loop110

  rax_av  =  0.0_DP
  zax_av  =  0.0_DP
  bax_av  =  0.0_DP
  bpax_av =  0.0_DP
  loop210 : DO i=1,lnx
    phi =  x0 + h * (i - 1)
    CALL mgval1(f(1,i), phi, f(2,i), br, bp, bz, bb)
    rax_av  =  rax_av  + f(1,i)
    zax_av  =  zax_av  + f(2,i)
    bax_av  =  bax_av  + bb
    bpax_av =  bpax_av + bp
  END DO loop210

  rax_av  =  rax_av  / lnx
  zax_av  =  zax_av  / lnx
  bax_av  =  bax_av  / lnx
  bpax_av =  bpax_av / lnx

  fmt =  '(//3(A7, F9.4), A8, F9.4//)'
  PRINT fmt, ' <R> = ', rax_av, ' <Z> = ', zax_av, ' <B> = ', bax_av, ' <Bp> = ', bpax_av

  DEALLOCATE(f0, f)


  RETURN
END SUBROUTINE mgaxis
