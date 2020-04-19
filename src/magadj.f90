SUBROUTINE magadj

  USE kind_spec
  USE axis_mod,              ONLY : rmagax,  &
    &                               fbvadj,  &
    &                               btor0,   &
    &                               rbtor0,  &
    &                               rax_av,  &
    &                               bpax_av
  USE cylindrical_coord_mod, ONLY : cj,      &
    &                               magset,  &
    &                               ipfcoil, &
    &                               cfact

  IMPLICIT NONE

  INTEGER :: i, &
    &        n
  REAL(DP) :: cpf_max, &
    &         cpf_min, &
    &         cpf_mid, &
    &         cpf_old


  rmagax =  3.60_DP
  fbvadj =  1.03_DP
  btor0  =  3.0_DP

  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE magadj                                      '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *

  cfact(:)   =  1.0_DP

  ipfcoil(:) =  0
  ipfcoil(4) =  1
  ipfcoil(5) =  1
  ipfcoil(6) =  1
  
  CALL mgaxis

  cpf_old =  1.0_DP

  DO i=1,50
    IF(ipfcoil(i) == 1)THEN
      cfact(i) =  fbvadj
    END IF
  END DO

  CALL magset

  CALL mgaxis

  IF(rax_av > rmagax)THEN
    cpf_min =  cpf_old
    cpf_max =  fbvadj
  ELSE
    cpf_min =  fbvadj
    cpf_max =  cpf_old
  END IF

  DO n=1,10!0

    PRINT *
    PRINT *, ' n= ', n
    PRINT *

    cpf_mid = (cpf_min + cpf_max) / 2

    DO i=1,50
      IF(ipfcoil(i) == 1)THEN
        cfact(i) =  cpf_mid
      END IF
    END DO

    CALL magset

    CALL mgaxis 

    IF(ABS(rax_av - rmagax) < 0.001_DP) EXIT

    IF(rax_av > rmagax)THEN
      cpf_max =  cpf_mid
    ELSE
      cpf_min =  cpf_mid
    END IF

  END DO

  CALL magset

  CALL mgaxis

  cfact(:) =  cfact(:) * btor0 / bpax_av

  PRINT *
  PRINT *, cj(1:6) * cfact(1:6)
  PRINT *

  CALL magset

  CALL mgaxis


  STOP  
END SUBROUTINE magadj
