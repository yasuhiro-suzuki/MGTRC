!> @file order.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: order
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
SUBROUTINE order (mplot,       & !(in)
  &               rval, zval,  & !(inout)
  &               raxis, zaxis & !(in)
  &              )

  USE kind_spec
  USE param1,   ONLY : pi2

  IMPLICIT NONE

!Arguments
  INTEGER, INTENT(IN) :: mplot
  REAL(DP), INTENT(IN) :: raxis, &
    &                     zaxis
  REAL(DP), INTENT(INOUT) :: rval(mplot), &
    &                        zval(mplot)
!Local variables
  INTEGER :: i,    &
    &        j,    &
    &        ip1,  &
    &        isave
  REAL(DP) :: ang1,  &
    &         ang2,  &
    &         r,     &
    &         z,     &
    &         saver, &
    &         savez


  loop010 : DO i=1,mplot-1
    ip1   =  i + 1
    isave =  i
    r     =  rval(i) - raxis
    z     =  zval(i) - zaxis
    ang1  =  ATAN2(z, r)
    IF(ang1 < 0.0_DP) ang1 =  ang1 + pi2
    IF(ang1 > pi2)    ang1 =  ang1 - pi2
    loop020 : DO j=ip1,mplot
      r    =  rval(j) - raxis
      z    =  zval(j) - zaxis
      ang2 =  ATAN2(z, r)
      IF(ang2 < 0.0_DP) ang2 =  ang2 + pi2
      IF(ang2 > pi2)    ang2 =  ang2 - pi2
      IF(ang2 < ang1)THEN
        isave =  j
        ang1  =  ang2
      END IF
    END DO loop020
    saver       =  rval(i) 
    savez       =  zval(i)
    rval(i)     =  rval(isave)
    zval(i)     =  zval(isave)
    rval(isave) =  saver
    zval(isave) =  savez
  END DO loop010

 
  RETURN
END SUBROUTINE order
