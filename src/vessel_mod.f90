!=vessel_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! SUBROUTINE for vacuum vessel
!
!==Reference
!
! None
!
!==Error Handlings
!
! None
!
!==Known Bugs
!
! None
!
!==Note
!
! None
!
!==TODO
!
! None
!

MODULE vessel_mod

  USE kind_spec
  USE param1,   ONLY : pi, &
    &                  pi2

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: lvessel,       &
    &        lcheck_vessel, &
    &        lvessel_vtk,   &
    &        lvessel_txt
  INTEGER :: mtor, &
    &        nu,   &
    &        nv
  CHARACTER(LEN=10) :: vessel_model = ''
  CHARACTER(LEN=300) :: vessel_file = 'vessel.dat'

  REAL(DP) :: pi2m, &
    &         dtor
  REAL(DP), ALLOCATABLE :: xx(:,:), &
    &                      yy(:,:), &
    &                      zz(:,:), &
    &                      rr(:,:)

  PUBLIC :: lvessel,         &
    &       lcheck_vessel,   &
    &       lvessel_vtk,     &
    &       vessel_model,    &
    &       vessel_file,     &
    &       nu,              &
    &       nv,              &
    &       read_vessel,     &
    &       free_mem_vessel, &
    &       vessel,          &
    &       check_vessel,    &
    &       vessel_loss

CONTAINS

  SUBROUTINE make_mem_vessel

    IMPLICIT NONE


    ALLOCATE(xx(nu,nv), yy(nu,nv), zz(nu,nv), rr(nu,nv))


  END SUBROUTINE make_mem_vessel

  SUBROUTINE free_mem_vessel

    IMPLICIT NONE


    DEALLOCATE(xx, yy, zz, rr)


  END SUBROUTINE free_mem_vessel

  SUBROUTINE read_vessel

    IMPLICIT NONE


    SELECT CASE(TRIM(vessel_model))
      CASE("2d", "2D")
        CALL read_vessel_2d
      CASE("3d", "3D")
        CALL read_vessel_3d
      CASE("3d_new", "3D_new")
        CALL read_vessel_3d_new
      CASE("old", "OLD")
        CALL read_vessel_old
      CASE("ana", "ANA")
        CALL read_vessel_ana
      CASE DEFAULT
        CALL read_vessel_3d
    END SELECT

    IF(lvessel_vtk) CALL write_vessel_vtk


  END SUBROUTINE read_vessel

  SUBROUTINE read_vessel_3d_new

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(DP) :: phi


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / nv

    nv   =  nv + 1

    CALL make_mem_vessel

    DO j=1,nv-1
      DO i=1,nu
        READ(30,*) rr(i,j), zz(i,j)
      END DO
    END DO

    DO i=1,nu
      rr(i,nv) =  rr(i,1)
      zz(i,nv) =  zz(i,1)
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  rr(i,j) * COS(phi)
        yy(i,j) =  rr(i,j) * SIN(phi)
      END DO
    END DO


  END SUBROUTINE read_vessel_3d_new

  SUBROUTINE read_vessel_3d

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(DP) :: phi


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO j=1,nv
      DO i=1,nu
        READ(30,*) rr(i,j), zz(i,j)
      END DO
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  rr(i,j) * COS(phi)
        yy(i,j) =  rr(i,j) * SIN(phi)
      END DO
    END DO


  END SUBROUTINE read_vessel_3d

  SUBROUTINE read_vessel_old

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j


    READ(30,*) mtor, nu, nv

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO j=1,nv
      DO i=1,nu
        READ(30,*) xx(i,j), yy(i,j), zz(i,j), rr(i,j)
      END DO
    END DO


  END SUBROUTINE read_vessel_old

  SUBROUTINE read_vessel_2d

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(DP) :: phi
    REAL(DP), ALLOCATABLE :: r1(:), &
      &                      z1(:)

    READ(30,*) nu

    ALLOCATE(r1(nu), z1(nu))

    mtor =  1
    nv   =  120

    pi2m =  pi2  / mtor
    dtor =  pi2m / (nv - 1)

    CALL make_mem_vessel

    DO i=1,nu
      READ(30,*) r1(i), z1(i)
    END DO

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        xx(i,j) =  r1(i) * COS(phi)
        yy(i,j) =  r1(i) * SIN(phi)
        zz(i,j) =  z1(i)
        rr(i,j) =  r1(i)
      END DO
    END DO

    DEALLOCATE(r1, z1)


  END SUBROUTINE read_vessel_2d

  SUBROUTINE read_vessel_ana

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j
    REAL(DP), PARAMETER :: r0 =  3.0_DP, &
      &                    a0 =  1.2_DP
    REAL(DP) :: theta,  &
      &         dtheta, &
      &         phi,    &
      &         r,      &
      &         z



    mtor   =  1
    nv     =  121
    nu     =  361

    pi2m   =  pi2  / mtor
    dtor   =  pi2m / (nv - 1)
    dtheta =  pi2 / (nu - 1)

    CALL make_mem_vessel

    DO j=1,nv
      phi =  dtor * (j - 1)
      DO i=1,nu
        theta   =  dtheta * (i - 1)
        r       =  r0 + a0 * COS(theta)
        z       =       a0 * SIN(theta)
        xx(i,j) =  r * COS(phi)
        yy(i,j) =  r * SIN(phi)
        zz(i,j) =  z
        rr(i,j) =  r
      END DO
    END DO


  END SUBROUTINE read_vessel_ana

  SUBROUTINE write_vessel_vtk

    IMPLICIT NONE
!Local variables
    INTEGER :: i1, &
      &        i2, &
      &        i3, &
      &        i,  &
      &        j


    WRITE(300,'(A)') '# vtk DataFile Version 3.0'
    WRITE(300,'(A)') 'Unstructured Grid'
    WRITE(300,'(A)') 'ASCII'
    WRITE(300,*)
    WRITE(300,'(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(300,'(A,I12,A)') 'POINTS ', nu * nv, ' double'

    DO j=1,nv
      DO i=1,nu
        WRITE(300,'(3(ES15.7))') xx(i,j), yy(i,j), zz(i,j)
      END DO
    END DO

    WRITE(300,*)
    WRITE(300,'(A,2I12)') 'CELLS', 2 * (nv - 1) * (nu - 1), 6 * (nv - 1) * (nu - 1) + 2 * (nv - 1) * (nu - 1)

    DO j=1,nv-1
      DO i=1,nu-1
        i1 =  i  + (j - 1) * nu 
        i2 =  i1 + 1
        i3 =  i  + j * nu
        WRITE(300,'(I3,3I12)') 3, i1 - 1, i2 - 1, i3 - 1
        i1 =  i  + (j - 1) * nu +1
        i2 =  i  + j * nu
        i3 =  i2 + 1 
        WRITE(300,'(I3,3I12)') 3, i1 - 1, i2 - 1, i3 - 1
      END DO
    END DO

    WRITE(300,*)
    WRITE(300,'(A,I12)') 'CELL_TYPES', 2 * (nv - 1) * (nu - 1)

    DO j=1,nv-1
      DO i=1,nu-1
        WRITE(300,'(I2)') 5
        WRITE(300,'(I2)') 5
      END DO
    END DO

    WRITE(300,'(A,I12)') 'POINT_DATA', nv * nu
    WRITE(300,'(A)') 'SCALARS scalars float 1'
    WRITE(300,'(A)') 'LOOKUP_TABLE default'

    DO j=1,nv
      DO i=1,nu
        WRITE(300,'(ES15.7)') 1.0_DP
      END DO
    END DO


  END SUBROUTINE write_vessel_vtk

  SUBROUTINE write_vessel_txt

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j


    DO j=1,nv-1
      DO i=1,nu-1
        WRITE(200,'(3ES15.7)') xx(i,j), yy(i,j), zz(i,j)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,'(3ES15.7)') xx(i,j+1), yy(i,j+1), zz(i,j+1)
        WRITE(200,'(3ES15.7)') xx(i,j), yy(i,j), zz(i,j)
        WRITE(200,*)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,'(3ES15.7)') xx(i,j+1), yy(i,j+1), zz(i,j+1)
        WRITE(200,'(3ES15.7)') xx(i+1,j+1), yy(i+1,j+1), zz(i+1,j+1)
        WRITE(200,'(3ES15.7)') xx(i+1,j), yy(i+1,j), zz(i+1,j)
        WRITE(200,*)
      END DO
      WRITE(200,*)
      WRITE(200,*)
    END DO


  END SUBROUTINE write_vessel_txt
 
  SUBROUTINE vessel (phi,   & !(in)
    &                r1, z1 & !(out)
    &               )


    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: phi
    REAL(DP), INTENT(OUT) :: r1(nu), &
      &                      z1(nu)
!Local variables
    INTEGER :: iphi
    REAL(DP) :: phi0,  &
      &         phi1,  &
      &         dtor1, &
      &         alpha


    iphi  =  phi / pi2m
    phi1  =  phi - pi2m * iphi
    IF(phi1 <  0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)   phi1 =  phi1 - pi2m

    iphi  =  phi1 / dtor
    phi0  =  dtor * iphi
    dtor1 =  phi1 - phi0

    iphi  =  iphi + 1

    alpha =  dtor1 / dtor

    r1(:) =  rr(:,iphi) + alpha * (rr(:,iphi+1) - rr(:,iphi))
    z1(:) =  zz(:,iphi) + alpha * (zz(:,iphi+1) - zz(:,iphi))


  END SUBROUTINE vessel

  SUBROUTINE check_vessel (r, z, phi,     & !(in)
    &                      iout           & !(out)
    &                     )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: r,  &
      &                     z,  &
      &                     phi
!Local variables
    INTEGER :: iphi, &
      &        i
    REAL(DP) :: phi0,   &
      &         phi1,   &
      &         dtor1,  &
      &         alpha,  &
      &         theta,  &
      &         theta0, &
      &         theta1, &
      &         theta2, &
      &         dtheta, &
      &         r1(nu), &
      &         z1(nu)


    iphi  =  phi / pi2m
    phi1  =  phi - pi2m * iphi
    IF(phi1 <  0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)   phi1 =  phi1 - pi2m

    iphi  =  phi1 / dtor
    phi0  =  dtor * iphi
    dtor1 =  phi1 - phi0

    iphi  =  iphi + 1

    alpha =  dtor1 / dtor

    r1(:) =  rr(:,iphi) + alpha * (rr(:,iphi+1) - rr(:,iphi))
    z1(:) =  zz(:,iphi) + alpha * (zz(:,iphi+1) - zz(:,iphi))

    theta  =  0.0_DP
    theta0 =  ATAN2(z1(1) - z, r1(1) - r)
    IF(theta0 < 0.0_DP) theta0 =  theta0 + pi2
    DO i=2,nu
      theta1 =  ATAN2(z1(i) - z, r1(i) - r)
      theta2 =  theta1
      IF(theta1 < 0.0_DP) theta1 =  theta1 + pi2
      dtheta =  theta1 - theta0
      IF(dtheta >  pi) dtheta =  dtheta - pi2
      IF(dtheta < -pi) dtheta =  dtheta + pi2
      theta  =  theta + dtheta
      theta0 =  theta1
    END DO
    theta =  theta / pi2
    IF(theta >= 0.0_DP)THEN
      theta =  theta + 1.0e-06_DP
    ELSE
      theta =  theta - 1.0e-06_DP
    END IF
    iout =  theta

 
  END SUBROUTINE check_vessel

  SUBROUTINE vessel_loss (r0, z0, p0, r1, z1, p1, &!(in)
    &                     point                   &!(out)
    &                    )

    IMPLICIT NONE

    REAL(DP) :: r0, z0, p0, r1, z1, p1

    INTEGER :: iphi0, iphi1, intersect, i, j
    REAL(DP) :: x0, y0, x1, y1, phi0, phi1, rt_t, rt_u, rt_v, save0(3), save1(3), point(3), origin(3), direction(3), rt_vert1(3), rt_vert2(3), rt_vert3(3)


    iphi0 =  p0 / pi2m
    iphi1 =  p1 / pi2m
    phi0  =  p0 - pi2m * iphi0
    phi1  =  p1 - pi2m * iphi1
    IF(phi0 <  0.0_DP) phi0 =  phi0 + pi2m
    IF(phi0 >= pi2m)   phi0 =  phi0 - pi2m
    IF(phi1 <  0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)   phi1 =  phi1 - pi2m

    save0(1) =  r0
    save0(2) =  z0
    save0(3) =  phi0

    save1(1) =  r1
    save1(2) =  z1
    save1(3) =  phi1

    IF(phi0 > phi1)THEN
      r0   =  save1(1)
      z0   =  save1(2)
      phi0 =  save1(3)
      r1   =  save0(1)
      z1   =  save0(2)
      phi1 =  save0(3)
    END IF

    iphi0  =  phi0 / dtor + 1
    iphi1  =  phi1 / dtor + 1

    x0 =  r0 * COS(phi0)
    y0 =  r0 * SIN(phi0)
    x1 =  r1 * COS(phi1)
    y1 =  r1 * SIN(phi1)

    origin(1)    =  x0
    origin(2)    =  y0
    origin(3)    =  z0
    direction(1) =  x1 - x0
    direction(2) =  y1 - y0
    direction(3) =  z1 - z0

    !DO j=1,nv-1
    !  DO i=1,nu-1
    !    rt_vert1(1) =  xx(i,j)
    !    rt_vert1(2) =  yy(i,j)
    !    rt_vert1(3) =  zz(i,j)
    !    rt_vert2(1) =  xx(i+1,j)
    !    rt_vert2(2) =  yy(i+1,j)
    !    rt_vert2(3) =  zz(i+1,j)
    !    rt_vert3(1) =  xx(i,j+1)
    !    rt_vert3(2) =  yy(i,j+1)
    !    rt_vert3(3) =  zz(i,j+1)
    !    CALL rt_intersec(origin(:), direction(:), rt_vert1(:), rt_vert2(:), rt_vert3(:), rt_t, rt_u, rt_v, intersect)
    !    IF(intersect == 1)THEN
    !      point(:) =  origin(:) + rt_t * direction(:)
    !      WRITE(80,'(10e19.7)') point(1), point(2), point(3)
    !    END IF
    !    rt_vert1(1) =  xx(i+1,j)
    !    rt_vert1(2) =  yy(i+1,j)
    !    rt_vert1(3) =  zz(i+1,j)
    !    rt_vert2(1) =  xx(i,j+1)
    !    rt_vert2(2) =  yy(i,j+1)
    !    rt_vert2(3) =  zz(i,j+1)
    !    rt_vert3(1) =  xx(i+1,j+1)
    !    rt_vert3(2) =  yy(i+1,j+1)
    !    rt_vert3(3) =  zz(i+1,j+1)
    !    CALL rt_intersec(origin(:), direction(:), rt_vert1(:), rt_vert2(:), rt_vert3(:), rt_t, rt_u, rt_v, intersect)
    !    IF(intersect == 1)THEN
    !      point(:) =  origin(:) + rt_t * direction(:)
    !      WRITE(80,'(10e19.7)') point(1), point(2), point(3)
    !    END IF
    !  END DO
    !END DO
    DO i=1,nu-1
      rt_vert1(1) =  xx(i,iphi0)
      rt_vert1(2) =  yy(i,iphi0)
      rt_vert1(3) =  zz(i,iphi0)
      rt_vert2(1) =  xx(i+1,iphi0)
      rt_vert2(2) =  yy(i+1,iphi0)
      rt_vert2(3) =  zz(i+1,iphi0)
      rt_vert3(1) =  xx(i,iphi0+1)
      rt_vert3(2) =  yy(i,iphi0+1)
      rt_vert3(3) =  zz(i,iphi0+1)
      CALL rt_intersec(origin(:), direction(:), rt_vert1(:), rt_vert2(:), rt_vert3(:), rt_t, rt_u, rt_v, intersect)
      IF(intersect == 1)THEN
        point(:) =  origin(:) + rt_t * direction(:)
        RETURN
      END IF
      rt_vert1(1) =  xx(i+1,iphi0)
      rt_vert1(2) =  yy(i+1,iphi0)
      rt_vert1(3) =  zz(i+1,iphi0)
      rt_vert2(1) =  xx(i,iphi0+1)
      rt_vert2(2) =  yy(i,iphi0+1)
      rt_vert2(3) =  zz(i,iphi0+1)
      rt_vert3(1) =  xx(i+1,iphi0+1)
      rt_vert3(2) =  yy(i+1,iphi0+1)
      rt_vert3(3) =  zz(i+1,iphi0+1)
      CALL rt_intersec(origin(:), direction(:), rt_vert1(:), rt_vert2(:), rt_vert3(:), rt_t, rt_u, rt_v, intersect)
      IF(intersect == 1)THEN
        point(:) =  origin(:) + rt_t * direction(:)
        RETURN
      END IF
    END DO


  END SUBROUTINE vessel_loss

  SUBROUTINE rt_intersec(origin, direction, vertex1, vertex2, vertex3, &!(in)
    &                    t, u, v, intersect                            &!(out)
    &                   )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: origin(3), direction(3), vertex1(3), vertex2(3), vertex3(3)
    REAL(DP), INTENT(OUT) :: t, u, v
    INTEGER, INTENT(OUT) ::  intersect
    ! 1: ray intersects triangle
    ! 0: ray dose not intersects triangle

    REAL(DP), PARAMETER :: epsilon =  1.0E-08_DP
    REAL(DP) :: det, inv_det, edge1(3), edge2(3), tvec(3), pvec(3), qvec(3)
    REAL(DP), EXTERNAL :: dot


    edge1(:) =  vertex2(:) - vertex1(:)
    edge2(:) =  vertex3(:) - vertex1(:)

    pvec(1)  =  direction(2) * edge2(3) - direction(3) * edge2(2)
    pvec(2)  =  direction(3) * edge2(1) - direction(1) * edge2(3)
    pvec(3)  =  direction(1) * edge2(2) - direction(2) * edge2(1)

    det      =  DOT_PRODUCT(edge1(:), pvec(:))
    IF(det < epsilon)THEN
      intersect =  0
    ELSE
      tvec(:)   =  origin(:) - vertex1(:)
      u         =  DOT_PRODUCT(tvec(:), pvec(:))
      IF((u < 0.0_DP) .OR. (u > det))THEN
        intersect =  0
      ELSE

        qvec(1) =  tvec(2) * edge1(3) - tvec(3) * edge1(2)
        qvec(2) =  tvec(3) * edge1(1) - tvec(1) * edge1(3)
        qvec(3) =  tvec(1) * edge1(2) - tvec(2) * edge1(1)

        v       =  DOT_PRODUCT(direction(:), qvec(:))
        IF((v < 0.0_DP) .OR. (u + v > det))THEN
          intersect =  0
        ELSE
          t         =  DOT_PRODUCT(edge2(:), qvec(:))
          IF((t < 0.0_DP) .OR. (det < t))THEN
            intersect =  0
          ELSE
            intersect =  1

            inv_det   =  1.0_DP / det
            t         =  t * inv_det
            u         =  u * inv_det
            v         =  v * inv_det
          END IF
        END IF
      END IF
    END IF


  END SUBROUTINE rt_intersec

END MODULE vessel_mod
