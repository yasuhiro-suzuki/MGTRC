!=subf3d.f90
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
SUBROUTINE subf3d ( t, x,    & !(in)
  &                 xp, iout & !(out)
  &               )

  USE kind_spec
  USE cylindrical_coord_mod, ONLY : mgval1

  IMPLICIT NONE

!Arguments
  INTEGER, INTENT(OUT) :: iout
  REAL(DP), INTENT(IN) :: t,   &
    &                     x(3)
  REAL(DP), INTENT(OUT) :: xp(3)
!Local variables
  REAL(DP) :: r,   &
    &         z,   &
    &         phi, &
    &         br,  &
    &         bp,  &
    &         bz,  &
    &         bb


  iout =  0
  r    =  x(1)
  z    =  x(2)
  phi  =  x(3)

!----------------------------------------
  CALL mgval1(r, phi, z, br, bp, bz, bb)
!----------------------------------------

  IF(bb == 0.0_DP)THEN
    xp   =  0.0_DP
    iout =  1
    RETURN
  END IF

  xp(1) =  br / bb
  xp(2) =  bz / bb
  xp(3) =  bp / (r * bb)


  RETURN
END SUBROUTINE subf3d
