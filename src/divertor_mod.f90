!> @file divertor_mod.f90
!------------------------------------------------------------------------------
!
! MODULE: divertor_mod
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
MODULE divertor_mod

  USE kind_spec
  USE param1,   ONLY : pi2

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: ldivertor, &
    &        ldivertor_vtk

  INTEGER :: mtor =  1
  REAL(DP) :: pi2m

  INTEGER :: num_nodes    =  14, &
    &        num_elements =  12
  INTEGER, ALLOCATABLE :: def_elements(:,:)
  REAL(DP), ALLOCATABLE :: x_elements(:), &
    &                      y_elements(:), &
    &                      z_elements(:)

  PUBLIC :: ldivertor, ldivertor_vtk,     &
    &       num_nodes, &
    &       num_elements, &
    &       def_elements, &
    &       x_elements, &
    &       y_elements, &
    &       z_elements, &
    &       read_divertor,     &
    &       free_mem_divertor, &
    &       divertor_loss

CONTAINS

  SUBROUTINE make_mem_divertor

    IMPLICIT NONE


    ALLOCATE(def_elements(num_elements,3), x_elements(num_nodes), y_elements(num_nodes), z_elements(num_nodes))


  END SUBROUTINE make_mem_divertor

  SUBROUTINE free_mem_divertor

    IMPLICIT NONE


    DEALLOCATE(def_elements, x_elements, y_elements, z_elements)


  END SUBROUTINE free_mem_divertor

  SUBROUTINE read_divertor

    IMPLICIT NONE


    !CALL make_mem_divertor
    !CALL read_nodes
    !CALL read_elements
    CALL read_divertor_w7code

    IF(ldivertor_vtk) CALL write_divertor_vtk

    pi2m =  pi2 / mtor


  END SUBROUTINE read_divertor

  SUBROUTINE read_nodes

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        i


    DO i=1,num_nodes
      READ(32,*) itemp, x_elements(i), y_elements(i), z_elements(i)
    END DO


  END SUBROUTINE read_nodes

  SUBROUTINE read_elements

    IMPLICIT NONE

    INTEGER :: i,       &
      &        itemp(2)


    DO i=1,num_elements
      READ(33,*) itemp(1), def_elements(i,1), def_elements(i,2), def_elements(i,3)
    END DO


  END SUBROUTINE read_elements

  SUBROUTINE read_divertor_w7code

    IMPLICIT NONE

     INTEGER :: nu, nv
!Local variables
    INTEGER :: i1, i2, i3, i, &
      &        j, &
      &        k
    REAL(DP) :: r, phi, z

    READ(31,*)
    READ(31,*) nv, nu, mtor
    PRINT *, nv, nu, mtor

    num_nodes    =  nv * nu
    num_elements =  2 * (nv - 1) * (nu - 1)

    CALL make_mem_divertor

    k =  0
    DO j=1,nv
      READ(31,*) phi
      phi =  phi * pi2 / 360.0_DP
      DO i=1,nu
        k =  k + 1
        READ(31,*) r, z
        x_elements(k) =  1.0E-02_DP * r * COS(phi)
        y_elements(k) =  1.0E-02_DP * r * SIN(phi)
        z_elements(k) =  1.0E-02_DP * z
        !WRITE(32,*) x_elements(k), y_elements(k), z_elements(k)
      END DO
    END DO

    k =  0
    DO j=1,nv-1
      DO i=1,nu-1
        i1 =  i  + (j - 1) * nu
        i2 =  i1 + 1
        i3 =  i  + j * nu
        k  =  k + 1
        def_elements(k,1) =  i1 - 1
        def_elements(k,2) =  i2 - 1
        def_elements(k,3) =  i3 - 1
        i1 =  i  + (j - 1) * nu + 1
        i2 =  i  + j * nu
        i3 =  i2 + 1
        k  =  k + 1
        def_elements(k,1) =  i1 - 1
        def_elements(k,2) =  i2 - 1
        def_elements(k,3) =  i3 - 1
      END DO
    END DO


  END SUBROUTINE read_divertor_w7code

  SUBROUTINE write_divertor_vtk

    IMPLICIT NONE
!Local variables
    INTEGER :: i1, &
      &        i2, &
      &        i3, &
      &        i


    DO i=1,num_elements
      i1 =  def_elements(i,1) + 1
      i2 =  def_elements(i,2) + 1
      i3 =  def_elements(i,3) + 1
      WRITE(32,*) i1, i2, i3
      WRITE(201,'(3ES15.7)') x_elements(i1), y_elements(i1), z_elements(i1) 
      WRITE(201,'(3ES15.7)') x_elements(i2), y_elements(i2), z_elements(i2) 
      WRITE(201,'(3ES15.7)') x_elements(i3), y_elements(i3), z_elements(i3) 
      WRITE(201,'(3ES15.7)') x_elements(i1), y_elements(i1), z_elements(i1) 
      WRITE(201,*)
      WRITE(201,*)
    END DO

    WRITE(301,'(A)') '# vtk DataFile Version 3.0'
    WRITE(301,'(A)') 'Unstructured Grid'
    WRITE(301,'(A)') 'ASCII'
    WRITE(301,*)
    WRITE(301,'(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(301,'(A,I12,A)') 'POINTS ', num_nodes, ' double'

    DO i=1,num_nodes
      WRITE(301,'(3(1PE15.7))') x_elements(i), y_elements(i), z_elements(i)
    END DO

    WRITE(301,*)
    WRITE(301,'(A,2I12)') 'CELLS', num_elements, 4 * num_elements

    DO i=1,num_elements
     WRITE(301,'(I3,3I12)') 3, def_elements(i,1), def_elements(i,2), def_elements(i,3)
    END DO

    WRITE(301,*)
    WRITE(301,'(A,I12)') 'CELL_TYPES', num_elements

    DO i=1,num_elements
      WRITE(301,'(I2)') 5
    END DO

    WRITE(301,'(A,I12)') 'POINT_DATA', num_nodes
    WRITE(301,'(A)') 'SCALARS scalars float 1'
    WRITE(301,'(A)') 'LOOKUP_TABLE default'

    DO i=1,num_nodes
      WRITE(301,'(1PE15.7)') 5.0_DP
    END DO


  END SUBROUTINE write_divertor_vtk

  SUBROUTINE divertor_loss (r0, phi0, z0, r1, phi1, z1, &!(in)
    &                       intersect, point            &!(out)
    &                      )

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: intersect
    REAL(DP), INTENT(INOUT) :: r0,   &
      &                        r1,   &
      &                        phi0, &
      &                        phi1, &
      &                        z0,   &
      &                        z1
    REAL(DP), INTENT(OUT) :: point(3)

    INTEGER :: iphi0, iphi1, i1, i2, i3, i
    REAL(DP) :: x0, y0, x1, y1, lr, rt_t, rt_u, rt_v, origin(3), direction(3), rt_vert1(3), rt_vert2(3), rt_vert3(3)


    iphi0 =  phi0 / pi2m
    iphi1 =  phi1 / pi2m
    phi0  =  phi0 - pi2m * iphi0
    phi1  =  phi1 - pi2m * iphi1
    IF(phi0 <  0.0_DP) phi0 =  phi0 + pi2m
    IF(phi0 >= pi2m)   phi0 =  phi0 - pi2m
    IF(phi1 <  0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)   phi1 =  phi1 - pi2m

    x0 =  r0 * COS(phi0)
    y0 =  r0 * SIN(phi0)
    x1 =  r1 * COS(phi1)
    y1 =  r1 * SIN(phi1)

    origin(1)    =  x0
    origin(2)    =  y0
    origin(3)    =  z0
    lr           =  SQRT((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)
    direction(1) = (x1 - x0) / lr
    direction(2) = (y1 - y0) / lr
    direction(3) = (z1 - z0) / lr

    DO i=1,num_elements
      i1 =  def_elements(i,1) + 1
      i2 =  def_elements(i,2) + 1
      i3 =  def_elements(i,3) + 1
      rt_vert1(1) =  x_elements(i1)
      rt_vert1(2) =  y_elements(i1)
      rt_vert1(3) =  z_elements(i1)
      rt_vert2(1) =  x_elements(i2)
      rt_vert2(2) =  y_elements(i2)
      rt_vert2(3) =  z_elements(i2)
      rt_vert3(1) =  x_elements(i3)
      rt_vert3(2) =  y_elements(i3)
      rt_vert3(3) =  z_elements(i3)
      CALL rt_intersec(origin(:), direction(:), rt_vert1(:), rt_vert2(:), rt_vert3(:), rt_t, rt_u, rt_v, intersect)
      IF(intersect == 1)THEN
        point(:) =  origin(:) + rt_t * direction(:)
        WRITE(180,'(10E19.7)') point(1), point(2), point(3)
        RETURN
      END IF
    END DO


  END SUBROUTINE divertor_loss

  SUBROUTINE rt_intersec(origin, direction, vertex1, vertex2, vertex3, &!(in)
    &                    t, u, v, intersect                            &!(out)
    &                   )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: origin(3), direction(3), vertex1(3), vertex2(3), vertex3(3)
    REAL(DP), INTENT(OUT) :: t, u, v
    INTEGER, INTENT(OUT) ::  intersect
    ! 1: ray intersects triangle
    ! 0: ray dose not intersects triangle

    REAL(DP), PARAMETER :: epsilon =  1.0E-01_DP
    REAL(DP) :: det, inv_det, edge1(3), edge2(3), tvec(3), pvec(3), qvec(3)
  

    edge1(:) =  vertex2(:) - vertex1(:)
    edge2(:) =  vertex3(:) - vertex1(:)

    pvec(1)  =  direction(2) * edge2(3) - direction(3) * edge2(2)
    pvec(2)  =  direction(3) * edge2(1) - direction(1) * edge2(3)
    pvec(3)  =  direction(1) * edge2(2) - direction(2) * edge2(1)

    det      =  DOT_PRODUCT(edge1(:), pvec(:))
    IF(det > epsilon)THEN
      tvec(:)   =  origin(:) - vertex1(:)
      u         =  DOT_PRODUCT(tvec(:), pvec(:))
      IF((u < 0.0_DP) .OR. (u > det))THEN
        intersect =  0
        RETURN
      END IF

      qvec(1) =  tvec(2) * edge1(3) - tvec(3) * edge1(2)
      qvec(2) =  tvec(3) * edge1(1) - tvec(1) * edge1(3)
      qvec(3) =  tvec(1) * edge1(2) - tvec(2) * edge1(1)

      v       =  DOT_PRODUCT(direction(:), qvec(:))
      IF((v < 0.0_DP) .OR. (u + v > det))THEN
        intersect =  0
        RETURN
      END IF
    ELSE IF(det < -epsilon)THEN
      tvec(:)   =  origin(:) - vertex1(:)
      u         =  DOT_PRODUCT(tvec(:), pvec(:))
      IF((u > 0.0_DP) .OR. (u < det))THEN
        intersect =  0
        RETURN
      END IF

      qvec(1) =  tvec(2) * edge1(3) - tvec(3) * edge1(2)
      qvec(2) =  tvec(3) * edge1(1) - tvec(1) * edge1(3)
      qvec(3) =  tvec(1) * edge1(2) - tvec(2) * edge1(1)

      v       =  DOT_PRODUCT(direction(:), qvec(:))
      IF((v > 0.0_DP) .OR. (u + v < det))THEN
        intersect =  0
        RETURN
      END IF

    ELSE
      intersect =  0
      RETURN
    END IF

    intersect =  1

    inv_det   =  1.0_DP / det
    t         =  t * inv_det
    u         =  u * inv_det
    v         =  v * inv_det
 

  END SUBROUTINE rt_intersec

END MODULE divertor_mod
