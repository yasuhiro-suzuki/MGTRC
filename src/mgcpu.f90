!> @file mgcpu.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: mgcpu
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
SUBROUTINE mgcpu (itime, mess & ! (in)
  &              )

  USE kind_spec

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: itime
  CHARACTER(LEN=*), INTENT(IN) :: mess

  INTEGER :: i,     &
    &        min1,  &
    &        min2
  INTEGER, SAVE :: it = 0
  REAL(DP) :: ccp,   &
    &         sec1,  &
    &         sec2,  &
    &         etime, &
    &         ta(2)
  REAL(DP), SAVE :: cpu1(100) =  0.0_DP, &
    &               cpu2(100) =  0.0_DP
  CHARACTER(LEN=20) :: mess1(100) = ''
  CHARACTER(LEN=100) :: fmt



  it  = it + 1
!  ccp = etime(ta)
  CALL CPU_TIME(ccp)

  IF(it == 1)THEN
    cpu1(it) =  ccp
    cpu2(it) =  ccp
  ELSE
    cpu1(it) =  ccp
    cpu2(it) =  ccp - cpu1(it-1)
  END IF

  mess1(it) =  mess

  IF(itime == 0) RETURN

  PRINT *

  fmt = '(I4, 2X, A20, A15, I5, A4, F9.3, A4, A11, I5, A4, F9.3, A4)'

  loop010 : DO  i=1,it
    min1 =  cpu1(i) / 60
    min2 =  cpu2(i) / 60
    sec1 =  cpu1(i) - min1 * 60
    sec2 =  cpu2(i) - min2 * 60
    PRINT fmt, i, mess1(i), '; Exec Time ...', min2, ' min', sec2, ' sec', &
      &                     '; Total ...',     min1, ' min', sec1, ' sec'
  END DO loop010


END SUBROUTINE mgcpu
