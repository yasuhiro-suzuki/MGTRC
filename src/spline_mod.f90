!=spline_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!    This subroutine is written by Yasuhiro Suzuki
!        at National Institute for Fusion Science ( NIFS )
!         2010/12/12
!
!    Based program is wrtten by T.Watanabe
!        at NIFS
!         1990/8/29
!
!  Multi-Dimensional Spline Interpolation  ( 4th )
!
!  This module contains following public routines.
!
!   1. splin1 : Initializing 1-dimensional interpolation.
!   2. splin2 : Initializing 2-dimensional interpolation.
!   3. splin3 : Initializing 3-dimensional interpolation.
!   4. sp1df  : 1-dimensional interpolation. 
!   5. sp1dd  : 1-dimensional interpolation and derivation along X. 
!   6. sp2df  : 2-dimensional interpolation. 
!   7. sp2dd  : 2-dimensional interpolation and derivation along X and Y.
!   8. sp3df  : 3-dimensional interpolation.
!   9. sp3dd  : 3-dimensional interpolation and derivation along X, Y and Z.
!
!
!              xsc                          xlc
!               <------------- x ------------>
!      1--+--+--4-------------------------nxxm-3-+--+--nxxm  a(i)
!      +--------+----------------------------+---------+
!               <---------------------------->
!
!
!
!==Reference
!
!  "Multi-Dimension Highly Accurate Spline Interpolation Method"
!  by T.Watanabe  The Japan Society for Industrial and Applied Mathmatics
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!
MODULE spline_mod

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: DP =  SELECTED_REAL_KIND(15)

  INTEGER :: l1d =  1, &
    &        l2d =  1, &
    &        l3d =  3, &
    &        nx1d,     &
    &        nx2d,     &
    &        ny2d,     &
    &        nx3d,     &
    &        ny3d,     &
    &        nz3d
   
  REAL(DP), PARAMETER :: c411 = -5.0_DP      / 2048.0_DP,   &
    &                    c412 =  9611.0_DP   / 737280.0_DP, &
    &                    c413 =  259.0_DP    / 23040.0_DP,  &
    &                    c414 = -6629.0_DP   / 46080.0_DP,  &
    &                    c415 = -7.0_DP      / 1152.0_DP,   &
    &                    c416 =  819.0_DP    / 1024.0_DP,   &
    &                    c417 =  1.0_DP      / 1440.0_DP,   &
    &                    c418 = -1067.0_DP   / 360.0_DP,    &
    &                    c41a =  3941.0_DP   / 576.0_DP,    &
    &                    c41c = -1603.0_DP   / 180.0_DP,    &
    &                    c41e =  901.0_DP    / 180.0_DP,    &
    &                    c421 =  49.0_DP     / 2048.0_DP,   &
    &                    c422 = -70733.0_DP  / 737280.0_DP, &
    &                    c423 = -499.0_DP    / 4608.0_DP,   &
    &                    c424 =  47363.0_DP  / 46080.0_DP,  &
    &                    c425 =  59.0_DP     / 1152.0_DP,   &
    &                    c426 = -86123.0_DP  / 15360.0_DP,  &
    &                    c427 = -1.0_DP      / 288.0_DP,    &
    &                    c431 = -245.0_DP    / 2048.0_DP,   &
    &                    c432 =  27759.0_DP  / 81920.0_DP,  &
    &                    c433 =  1299.0_DP   / 2560.0_DP,   &
    &                    c434 = -50563.0_DP  / 15360.0_DP,  &
    &                    c435 = -15.0_DP     / 128.0_DP,    &
    &                    c436 =  51725.0_DP  / 3072.0_DP,   &
    &                    c437 =  1.0_DP      / 160.0_DP,    &
    &                    c441 =  1225.0_DP   / 2048.0_DP,   &
    &                    c442 = -240077.0_DP / 147456.0_DP, &
    &                    c443 = -1891.0_DP   / 4608.0_DP,   &
    &                    c444 =  52931.0_DP  / 9216.0_DP,   &
    &                    c445 =  83.0_DP     / 1152.0_DP,   &
    &                    c446 = -86251.0_DP  / 3072.0_DP,   &
    &                    c447 = -1.0_DP      / 288.0_DP,    &
    &                    d413 =  c413 * 2.0_DP,                &
    &                    d423 =  c423 * 2.0_DP,                &
    &                    d433 =  c433 * 2.0_DP,                &
    &                    d443 =  c443 * 2.0_DP,                &
    &                    d414 =  c414 * 3.0_DP,                &
    &                    d424 =  c424 * 3.0_DP,                &
    &                    d434 =  c434 * 3.0_DP,                &
    &                    d444 =  c444 * 3.0_DP,                &
    &                    d415 =  c415 * 4.0_DP,                &
    &                    d425 =  c425 * 4.0_DP,                &
    &                    d435 =  c435 * 4.0_DP,                &
    &                    d445 =  c445 * 4.0_DP,                &
    &                    d416 =  c416 * 5.0_DP,                &
    &                    d426 =  c426 * 5.0_DP,                &
    &                    d436 =  c436 * 5.0_DP,                &
    &                    d446 =  c446 * 5.0_DP,                &
    &                    d417 =  c417 * 6.0_DP,                &
    &                    d427 =  c427 * 6.0_DP,                &
    &                    d437 =  c437 * 6.0_DP,                &
    &                    d447 =  c447 * 6.0_DP,                &
    &                    d418 =  c418 * 7.0_DP,                &
    &                    d41a =  c41a * 9.0_DP,                &
    &                    d41c =  c41c * 11.0_DP,               &
    &                    d41e =  c41e * 13.0_DP

  REAL(DP) :: h1x, &
    &         h2x, &
    &         h2y, &
    &         h3x, &
    &         h3y, &
    &         h3z, &
    &         xs1, &
    &         xs2, &
    &         ys2, &
    &         xs3, &
    &         ys3, &
    &         zs3, &
    &         xl1, &
    &         xl2, &
    &         yl2, &
    &         xl3, &
    &         yl3, &
    &         zl3
  REAL(DP), ALLOCATABLE :: f1d(:,:),    &
    &                      f2d(:,:,:),  &
    &                      f3d(:,:,:,:)
  !$acc declare create (f3d)

  PUBLIC :: l1d,    &
    &       l2d,    &
    &       l3d,    &
    &       nx1d,   &
    &       nx2d,   &
    &       ny2d,   &
    &       nx3d,   &
    &       ny3d,   &
    &       nz3d,   &
    &       f1d,    &
    &       f2d,    &
    &       f3d,    &
    &       splin1, &
    &       splin2, &
    &       splin3, &
    &       spl1df, &
    &       spl1dd, &
    &       spl2df, &
    &       spl2dd, &
    &       spl3df, &
    &       spl3dd

CONTAINS

  SUBROUTINE splin1 (xsd, xld & ! (in)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld


    xs1 =  xsd
    xl1 =  xld
    h1x = (xld - xsd) / (nx1d - 1)


  END SUBROUTINE splin1

  SUBROUTINE splin2 (xsd, xld, ysd, yld & ! (in)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld, &
      &                     ysd, &
      &                     yld


    xs2 =  xsd
    xl2 =  xld
    h2x = (xld - xsd) / (nx2d - 1)
    ys2 =  ysd
    yl2 =  yld
    h2y = (yld - ysd) / (ny2d - 1)


  END SUBROUTINE splin2

  SUBROUTINE splin3 (xsd, xld, ysd, yld, zsd, zld & ! (in)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld, &
      &                     ysd, &
      &                     yld, &
      &                     zsd, &
      &                     zld


    xs3 =  xsd  
    xl3 =  xld  
    h3x = (xld - xsd) / (nx3d - 1)
    ys3 =  ysd  
    yl3 =  yld 
    h3y = (yld - ysd) / (ny3d - 1)
    zs3 =  zsd  
    zl3 =  zld  
    h3z = (zld - zsd) / (nz3d - 1)


  END SUBROUTINE splin3

  SUBROUTINE spl1df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd
    REAL(DP), INTENT(OUT) :: w0(l1d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l, &
      &        ix
    REAL(DP) :: x,    &
      &         ux,   &
      &         usx,  &
      &         x400, &
      &         x41m, &
      &         x41p, &
      &         x42m, &
      &         x42p, &
      &         x43m, &
      &         x43p, &
      &         x44m, &
      &         x44p


    w0(:) =  0.0_DP


    x =  xd
    IF(x < xs1) RETURN
    IF(x > xl1) RETURN

    ux = (x - xs1) / h1x
    ix =  ux
    IF(ix >= nx1d - 1)THEN
      ix =  nx1d - 2
      ux =  nx1d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    usx   =  ux * ux

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx

    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441

    loop010 : DO l=1,l1d
      w0(l) = ((x41p + x41m) * f1d(l,ix-2)   &
        &   +  (x42p + x42m) * f1d(l,ix-1)   &
        &   +  (x43p + x43m) * f1d(l,ix  )   &
        &   +  (x44p + x44m) * f1d(l,ix+1)   &
        &   +  (x44p - x44m) * f1d(l,ix+2)   &
        &   +  (x43p - x43m) * f1d(l,ix+3)   &
        &   +  (x42p - x42m) * f1d(l,ix+4)   &
        &   +  (x41p - x41m) * f1d(l,ix+5))
    END DO loop010


  END SUBROUTINE spl1df

  SUBROUTINE spl1dd (xd,    &!(in)
    &                w0, wx &!(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l1d), &
      &                      wx(l1d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l, &
      &        ix
    REAL(DP) :: x,     &
      &         ux,    &
      &         usx,   &
      &         x400,  &
      &         x41m,  &
      &         x41p,  &
      &         x42m,  &
      &         x42p,  &
      &         x43m,  &
      &         x43p,  &
      &         x44m,  &
      &         x44p,  &
      &         dx400, &
      &         dx41m, &
      &         dx41p, &
      &         dx42m, &
      &         dx42p, &
      &         dx43m, &
      &         dx43p, &
      &         dx44m, &
      &         dx44p


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs1) RETURN
    IF(x > xl1) RETURN

    ux = (x - xs1) / h1x
    ix =  ux
    IF(ix >= nx1d - 1)THEN
      ix =  nx1d - 2
      ux =  nx1d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    usx   =  ux * ux

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux

    x42m  = (((-x400  * 7.0_DP  + c426) * usx + c424) * usx + c422) * ux
    x42p  =  ((c427   * usx     + c425) * usx + c423) * usx + c421
    dx42m =  ((-dx400 * 7.0_DP  + d426) * usx + d424) * usx + c422
    dx42p =  ((d427   * usx     + d425) * usx + d423) * ux

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux

    loop010 : DO l=1,l1d
      w0(l) = ((x41p + x41m) * f1d(l,ix-2)   &
        & +    (x42p + x42m) * f1d(l,ix-1)   &
        & +    (x43p + x43m) * f1d(l,ix  )   &
        & +    (x44p + x44m) * f1d(l,ix+1)   &
        & +    (x44p - x44m) * f1d(l,ix+2)   &
        & +    (x43p - x43m) * f1d(l,ix+3)   &
        & +    (x42p - x42m) * f1d(l,ix+4)   &
        & +    (x41p - x41m) * f1d(l,ix+5))

      wx(l) = ((dx41p + dx41m) * f1d(l,ix-2)  &
        & +    (dx42p + dx42m) * f1d(l,ix-1)  &
        & +    (dx43p + dx43m) * f1d(l,ix  )  &
        & +    (dx44p + dx44m) * f1d(l,ix+1)  &
        & +    (dx44p - dx44m) * f1d(l,ix+2)  &
        & +    (dx43p - dx43m) * f1d(l,ix+3)  &
        & +    (dx42p - dx42m) * f1d(l,ix+4)  &
        & +    (dx41p - dx41m) * f1d(l,ix+5))

      wx(l) =  wx(l) / h1x

    END DO loop010


  END SUBROUTINE spl1dd

  SUBROUTINE spl2df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l2d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l,  &
      &        my, &
      &        ly, &
      &        ix, &
      &        iy
    REAL(DP) :: x,     &
      &         y,     &
      &         ux,    &
      &         uy,    &
      &         usx,   &
      &         usy,   &
      &         x400,  &
      &         y400,  &
      &         x41m,  &
      &         y41m,  &
      &         x41p,  &
      &         y41p,  &
      &         x42m,  &
      &         y42m,  &
      &         x42p,  &
      &         y42p,  &
      &         x43m,  &
      &         y43m,  &
      &         x43p,  &
      &         y43p,  &
      &         x44m,  &
      &         y44m,  &
      &         x44p,  &
      &         y44p,  &
      &         cy(8)


    w0(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs2) RETURN
    IF(x > xl2) RETURN

    y =  xd(2)   
    IF(y < ys2) RETURN
    IF(y > yl2) RETURN

    ux = (x - xs2) / h2x
    ix =  ux
    IF(ix >= nx2d - 1)THEN
      ix =  nx2d - 2
      ux =  nx2d - 1
    END IF
 
    uy = (y - ys2) / h2y
    iy =  uy
    IF(iy >= ny2d - 1)THEN
      iy =  ny2d - 2
      uy =  ny2d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    uy    =  uy - iy - 0.5_DP
    usx   =  ux * ux
    usy   =  uy * uy

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400 = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    
    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41m = (((y400       + c416) * usy + c414) * usy + c412) * uy
    y41p =  ((c417 * usy + c415) * usy + c413) * usy + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421
    y42m = (((-y400 * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    y42p =  ((c427  * usy    + c425) * usy + c423) * usy + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431
    y43m = (((y400 * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    y43p =  ((c437 * usy     + c435) * usy + c433) * usy + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441
    y44m = (((-y400 * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    y44p =  ((c447  * usy     + c445) * usy + c443) * usy + c441

    cy(1) =  y41p + y41m
    cy(2) =  y42p + y42m
    cy(3) =  y43p + y43m
    cy(4) =  y44p + y44m
    cy(5) =  y44p - y44m
    cy(6) =  y43p - y43m
    cy(7) =  y42p - y42m
    cy(8) =  y41p - y41m

    loop010 : DO l=1,l2d
      loop020 : DO my=1,m2
        ly    =  iy + my - 3
        w0(l) = ((x41p + x41m) * f2d(l,ix-2,ly)                    &
          &   +  (x42p + x42m) * f2d(l,ix-1,ly)                    &
          &   +  (x43p + x43m) * f2d(l,ix  ,ly)                    &
          &   +  (x44p + x44m) * f2d(l,ix+1,ly)                    &
          &   +  (x44p - x44m) * f2d(l,ix+2,ly)                    & 
          &   +  (x43p - x43m) * f2d(l,ix+3,ly)                    &
          &   +  (x42p - x42m) * f2d(l,ix+4,ly)                    &
          &   +  (x41p - x41m) * f2d(l,ix+5,ly)) * cy(my) + w0(l)
      END DO loop020
    END DO loop010


  END SUBROUTINE spl2df

  SUBROUTINE spl2dd (xd,        &!(in)
    &                w0, wx, wy &!(out)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l2d), &
      &                      wx(l2d), &
      &                      wy(l2d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l,  &
      &        ix, &
      &        iy, &
      &        my, &
      &        ly
    REAL(DP) :: x,     &
      &         y,     &
      &         ux,    &
      &         uy,    &
      &         usx,   &
      &         usy,   &
      &         x400,  &
      &         y400,  &
      &         x41m,  &
      &         y41m,  &
      &         x41p,  &
      &         y41p,  &
      &         x42m,  &
      &         y42m,  &
      &         x42p,  &
      &         y42p,  &
      &         x43m,  &
      &         y43m,  &
      &         x43p,  &
      &         y43p,  &
      &         x44m,  &
      &         y44m,  &
      &         x44p,  &
      &         y44p,  &
      &         dx400, &
      &         dy400, &
      &         dx41m, &
      &         dy41m, &
      &         dx41p, &
      &         dy41p, &
      &         dx42m, &
      &         dy42m, &
      &         dx42p, &
      &         dy42p, &
      &         dx43m, &
      &         dy43m, &
      &         dx43p, &
      &         dy43p, &
      &         dx44m, &
      &         dy44m, &
      &         dx44p, &
      &         dy44p, &
      &         w0l,   &
      &         wxl,   &
      &         wyl,   &
      &         ww,    &
      &         cy(8), &
      &         dy(8)


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP
    wy(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs2) RETURN
    IF(x > xl2) RETURN
    
    y =  xd(2)
    IF(y < ys2) RETURN
    IF(y > yl2) RETURN

    ux = (x - xs2) / h2x
    ix =  ux
    IF(ix >= nx2d - 1)THEN
      ix =  nx2d - 2
      ux =  nx2d - 1
    END IF
 
    uy = (y - ys2) / h2y
    iy =  uy
    IF(iy >= ny2d - 1)THEN
      iy =  ny2d - 2
      uy =  ny2d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP

    usx =  ux * ux
    usy =  uy * uy

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400  = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx
    dy400 = (((d41e * usy + d41c) * usy + d41a) * usy + d418) * usy

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    y41m  = (((y400       + c416) * usy + c414) * usy + c412) * uy
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41p  =  ((c417 * usy + c415) * usy + c413) * usy + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dy41m =  ((dy400      + d416) * usy + d414) * usy + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux
    dy41p =  ((d417 * usy + d415) * usy + d413) * uy

    x42m  = (((-x400  * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    y42m  = (((-y400  * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    x42p  =  ((c427   * usx    + c425) * usx + c423) * usx + c421
    y42p  =  ((c427   * usy    + c425) * usy + c423) * usy + c421
    dx42m =  ((-dx400 * 7.0_DP + d426) * usx + d424) * usx + c422
    dy42m =  ((-dy400 * 7.0_DP + d426) * usy + d424) * usy + c422
    dx42p =  ((d427   * usx    + d425) * usx + d423) * ux
    dy42p =  ((d427   * usy    + d425) * usy + d423) * uy

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    y43m  = (((y400  * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    y43p  =  ((c437  * usy     + c435) * usy + c433) * usy + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dy43m =  ((dy400 * 21.0_DP + d436) * usy + d434) * usy + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux
    dy43p =  ((d437  * usy     + d435) * usy + d433) * uy

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    y44m  = (((-y400  * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    y44p  =  ((c447   * usy     + c445) * usy + c443) * usy + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dy44m =  ((-dy400 * 35.0_DP + d446) * usy + d444) * usy + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux
    dy44p =  (( d447  * usy     + d445) * usy + d443) * uy

    cy(1)   =  y41p +  y41m
    cy(2)   =  y42p +  y42m
    cy(3)   =  y43p +  y43m
    cy(4)   =  y44p +  y44m
    cy(5)   =  y44p -  y44m
    cy(6)   =  y43p -  y43m
    cy(7)   =  y42p -  y42m
    cy(8)   =  y41p -  y41m
    dy(1)   =  dy41p + dy41m
    dy(2)   =  dy42p + dy42m
    dy(3)   =  dy43p + dy43m
    dy(4)   =  dy44p + dy44m
    dy(5)   =  dy44p - dy44m
    dy(6)   =  dy43p - dy43m
    dy(7)   =  dy42p - dy42m
    dy(8)   =  dy41p - dy41m


    loop010 : DO l=1,l2d
      w0l =  0.0_DP
      wxl =  0.0_DP
      wyl =  0.0_DP

      loop020 : DO my=1,m2
        ly  =  iy + my - 3
        ww  = ((x41p + x41m) * f2d(l,ix-2,ly)   &
          & +  (x42p + x42m) * f2d(l,ix-1,ly)   &
          & +  (x43p + x43m) * f2d(l,ix  ,ly)   &
          & +  (x44p + x44m) * f2d(l,ix+1,ly)   &
          & +  (x44p - x44m) * f2d(l,ix+2,ly)   &
          & +  (x43p - x43m) * f2d(l,ix+3,ly)   &
          & +  (x42p - x42m) * f2d(l,ix+4,ly)   &
          & +  (x41p - x41m) * f2d(l,ix+5,ly))

        w0l =  ww * cy(my) + w0l
        wyl =  ww * dy(my) + wyl
        wxl = ((dx41p + dx41m) * f2d(l,ix-2,ly)                  &
          & +  (dx42p + dx42m) * f2d(l,ix-1,ly)                  &
          & +  (dx43p + dx43m) * f2d(l,ix  ,ly)                  &
          & +  (dx44p + dx44m) * f2d(l,ix+1,ly)                  &
          & +  (dx44p - dx44m) * f2d(l,ix+2,ly)                  &
          & +  (dx43p - dx43m) * f2d(l,ix+3,ly)                  &
          & +  (dx42p - dx42m) * f2d(l,ix+4,ly)                  &
          & +  (dx41p - dx41m) * f2d(l,ix+5,ly)) * cy(my) + wxl

      END DO loop020
         
      w0(l) =  w0l
      wx(l) =  wxl / h2x
      wy(l) =  wyl / h2y

    END DO loop010


    RETURN
  END SUBROUTINE spl2dd

  SUBROUTINE spl3df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd(3)
    REAL(DP), INTENT(OUT) :: w0(l3d)
!Local variables
    INTEGER, PARAMETER :: m2 = 8
    INTEGER :: l,  &
      &        mz, &
      &        lz, &
      &        my, &
      &        ly, &
      &        ix, &
      &        iy, &
      &        iz
    REAL(DP) :: x,     &
      &         y,     &
      &         z,     &
      &         ux,    &
      &         uy,    &
      &         uz,    &
      &         usx,   &
      &         usy,   &
      &         usz,   &
      &         x400,  &
      &         y400,  &
      &         z400,  &
      &         x41m,  &
      &         y41m,  &
      &         z41m,  &
      &         x41p,  &
      &         y41p,  &
      &         z41p,  &
      &         x42m,  &
      &         y42m,  &
      &         z42m,  &
      &         x42p,  &
      &         y42p,  &
      &         z42p,  &
      &         x43m,  &
      &         y43m,  &
      &         z43m,  &
      &         x43p,  &
      &         y43p,  &
      &         z43p,  &
      &         x44m,  &
      &         y44m,  &
      &         z44m,  &
      &         x44p,  &
      &         y44p,  &
      &         z44p,  &
      &         cyz,   &
      &         cy(8), &
      &         cz(8)


    w0(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs3) RETURN
    IF(x > xl3) RETURN

    y =  xd(2)   
    IF(y < ys3) RETURN
    IF(y > yl3) RETURN

    z =  xd(3)
    IF(z < zs3) z =  zs3
    IF(z > zl3) z =  zl3

    ux = (x - xs3) / h3x
    ix =  ux
    IF(ix >= nx3d - 1)THEN
      ix =  nx3d - 2
      ux =  nx3d - 1
    END IF
 
    uy = (y - ys3) / h3y
    iy =  uy
    IF(iy >= ny3d - 1)THEN
      iy =  ny3d - 2
      uy =  ny3d - 1
    END IF

    uz = (z - zs3) / h3z
    iz =  uz
    IF(iz >= nz3d - 1)THEN
      iz =  nz3d - 2
      uz =  nz3d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP
    uz  =  uz - iz - 0.5_DP
    usx =  ux * ux
    usy =  uy * uy
    usz =  uz * uz

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400 = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    z400 = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
    
    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41m = (((y400       + c416) * usy + c414) * usy + c412) * uy
    y41p =  ((c417 * usy + c415) * usy + c413) * usy + c411
    z41m = (((z400       + c416) * usz + c414) * usz + c412) * uz
    z41p =  ((c417 * usz + c415) * usz + c413) * usz + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421
    y42m = (((-y400 * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    y42p =  ((c427  * usy    + c425) * usy + c423) * usy + c421
    z42m = (((-z400 * 7.0_DP + c426) * usz + c424) * usz + c422) * uz
    z42p =  ((c427  * usz    + c425) * usz + c423) * usz + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431
    y43m = (((y400 * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    y43p =  ((c437 * usy     + c435) * usy + c433) * usy + c431
    z43m = (((z400 * 21.0_DP + c436) * usz + c434) * usz + c432) * uz
    z43p =  ((c437 * usz     + c435) * usz + c433) * usz + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441
    y44m = (((-y400 * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    y44p =  ((c447  * usy     + c445) * usy + c443) * usy + c441
    z44m = (((-z400 * 35.0_DP + c446) * usz + c444) * usz + c442) * uz
    z44p =  ((c447  * usz     + c445) * usz + c443) * usz + c441

    cy(1) =  y41p + y41m
    cy(2) =  y42p + y42m
    cy(3) =  y43p + y43m
    cy(4) =  y44p + y44m
    cy(5) =  y44p - y44m
    cy(6) =  y43p - y43m
    cy(7) =  y42p - y42m
    cy(8) =  y41p - y41m
    cz(1) =  z41p + z41m
    cz(2) =  z42p + z42m
    cz(3) =  z43p + z43m
    cz(4) =  z44p + z44m
    cz(5) =  z44p - z44m
    cz(6) =  z43p - z43m
    cz(7) =  z42p - z42m
    cz(8) =  z41p - z41m

    loop010 : DO l=1,l3d
      loop020 : DO mz=1,m2
        lz =  iz + mz - 3
        loop030 : DO my=1,m2
          ly    =  iy + my - 3
          cyz   =  cy(my) * cz(mz)
          w0(l) = ((x41p + x41m) * f3d(l,ix-2,ly,lz)                &
            &   +  (x42p + x42m) * f3d(l,ix-1,ly,lz)                &
            &   +  (x43p + x43m) * f3d(l,ix  ,ly,lz)                &
            &   +  (x44p + x44m) * f3d(l,ix+1,ly,lz)                &
            &   +  (x44p - x44m) * f3d(l,ix+2,ly,lz)                & 
            &   +  (x43p - x43m) * f3d(l,ix+3,ly,lz)                &
            &   +  (x42p - x42m) * f3d(l,ix+4,ly,lz)                &
            &   +  (x41p - x41m) * f3d(l,ix+5,ly,lz)) * cyz + w0(l)
        END DO loop030
      END DO loop020
    END DO loop010


    RETURN
  END SUBROUTINE spl3df

  SUBROUTINE spl3dd (xd,            &!(in)
    &                w0, wx, wy, wz &!(out)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(3)
    REAL(DP), INTENT(OUT) :: w0(l3d), &
      &                      wx(l3d), &
      &                      wy(l3d), &
      &                      wz(l3d)
!Local variables
    INTEGER, PARAMETER :: m2 = 8
    INTEGER :: l,   &
      &        ix,  &
      &        iy,  &
      &        iz,  &
      &        mz,  &
      &        lz,  &
      &        ix1, &
      &        ix2, &
      &        ix4, &
      &        ix5, &
      &        ix6, &
      &        ix7, &
      &        ix8, &
      &        iy1, &
      &        iy2, &
      &        iy4, &
      &        iy5, &
      &        iy6, &
      &        iy7, &
      &        iy8
    REAL(DP) :: x,     &
      &         y,     &
      &         z,     &
      &         ux,    &
      &         uy,    &
      &         uz,    &
      &         usx,   &
      &         usy,   &
      &         usz,   &
      &         x400,  &
      &         y400,  &
      &         z400,  &
      &         x41m,  &
      &         y41m,  &
      &         z41m,  &
      &         x41p,  &
      &         y41p,  &
      &         z41p,  &
      &         x42m,  &
      &         y42m,  &
      &         z42m,  &
      &         x42p,  &
      &         y42p,  &
      &         z42p,  &
      &         x43m,  &
      &         y43m,  &
      &         z43m,  &
      &         x43p,  &
      &         y43p,  &
      &         z43p,  &
      &         x44m,  &
      &         y44m,  &
      &         z44m,  &
      &         x44p,  &
      &         y44p,  &
      &         z44p,  &
      &         dx400, &
      &         dy400, &
      &         dz400, &
      &         dx41m, &
      &         dy41m, &
      &         dz41m, &
      &         dx41p, &
      &         dy41p, &
      &         dz41p, &
      &         dx42m, &
      &         dy42m, &
      &         dz42m, &
      &         dx42p, &
      &         dy42p, &
      &         dz42p, &
      &         dx43m, &
      &         dy43m, &
      &         dz43m, &
      &         dx43p, &
      &         dy43p, &
      &         dz43p, &
      &         dx44m, &
      &         dy44m, &
      &         dz44m, &
      &         dx44p, &
      &         dy44p, &
      &         dz44p, &
      &         cx1,   &
      &         cx2,   &
      &         cx3,   &
      &         cx4,   &
      &         cx5,   &
      &         cx6,   &
      &         cx7,   &
      &         cx8,   &
      &         cy1,   &
      &         cy2,   &
      &         cy3,   &
      &         cy4,   &
      &         cy5,   &
      &         cy6,   &
      &         cy7,   &
      &         cy8,   &
      &         dx1,   &
      &         dx2,   &
      &         dx3,   &
      &         dx4,   &
      &         dx5,   &
      &         dx6,   &
      &         dx7,   &
      &         dx8,   &
      &         dy1,   &
      &         dy2,   &
      &         dy3,   &
      &         dy4,   &
      &         dy5,   &
      &         dy6,   &
      &         dy7,   &
      &         dy8,   &
      &         w0l,   &
      &         wxl,   &
      &         wyl,   &
      &         wzl,   &
      &         w01,   &
      &         w02,   &
      &         w03,   &
      &         w04,   &
      &         w05,   &
      &         w06,   &
      &         w07,   &
      &         w08,   &
      &         w00,   &
      &         cz(8), &
      &         dz(8)


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP
    wy(:) =  0.0_DP
    wz(:) =  0.0_DP

    x =  xd(1)
    IF(x < xs3) RETURN
    IF(x > xl3) RETURN

    y =  xd(2)
    IF(y < ys3) RETURN
    IF(y > yl3) RETURN

    z =  xd(3)
    IF(z < zs3) z =  zs3
    IF(z > zl3) z =  zl3

    ux = (x - xs3) / h3x
    ix =  ux
    IF(ix >= nx3d - 1)THEN
      ix =  nx3d - 2
      ux =  nx3d - 1
    END IF
 
    uy = (y - ys3) / h3y
    iy =  uy
    IF(iy >= ny3d - 1)THEN
      iy =  ny3d - 2
      uy =  ny3d - 1
    END IF

    uz = (z - zs3) / h3z
    iz =  uz
    IF(iz >= nz3d - 1)THEN
      iz =  nz3d - 2
      uz =  nz3d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP
    uz  =  uz - iz - 0.5_DP

    usx =  ux * ux
    usy =  uy * uy
    usz =  uz * uz

    ix1 =  ix - 2
    ix2 =  ix - 1
    ix4 =  ix + 1
    ix5 =  ix + 2
    ix6 =  ix + 3
    ix7 =  ix + 4
    ix8 =  ix + 5

    iy1 =  iy - 2
    iy2 =  iy - 1
    iy4 =  iy + 1
    iy5 =  iy + 2
    iy6 =  iy + 3
    iy7 =  iy + 4
    iy8 =  iy + 5

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400  = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    z400  = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx
    dy400 = (((d41e * usy + d41c) * usy + d41a) * usy + d418) * usy
    dz400 = (((d41e * usz + d41c) * usz + d41a) * usz + d418) * usz

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    y41m  = (((y400       + c416) * usy + c414) * usy + c412) * uy
    z41m  = (((z400       + c416) * usz + c414) * usz + c412) * uz
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41p  =  ((c417 * usy + c415) * usy + c413) * usy + c411
    z41p  =  ((c417 * usz + c415) * usz + c413) * usz + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dy41m =  ((dy400      + d416) * usy + d414) * usy + c412
    dz41m =  ((dz400      + d416) * usz + d414) * usz + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux
    dy41p =  ((d417 * usy + d415) * usy + d413) * uy
    dz41p =  ((d417 * usz + d415) * usz + d413) * uz

    x42m  = (((-x400  * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    y42m  = (((-y400  * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    z42m  = (((-z400  * 7.0_DP + c426) * usz + c424) * usz + c422) * uz
    x42p  =  ((c427   * usx    + c425) * usx + c423) * usx + c421
    y42p  =  ((c427   * usy    + c425) * usy + c423) * usy + c421
    z42p  =  ((c427   * usz    + c425) * usz + c423) * usz + c421
    dx42m =  ((-dx400 * 7.0_DP + d426) * usx + d424) * usx + c422
    dy42m =  ((-dy400 * 7.0_DP + d426) * usy + d424) * usy + c422
    dz42m =  ((-dz400 * 7.0_DP + d426) * usz + d424) * usz + c422
    dx42p =  ((d427   * usx    + d425) * usx + d423) * ux
    dy42p =  ((d427   * usy    + d425) * usy + d423) * uy
    dz42p =  ((d427   * usz    + d425) * usz + d423) * uz

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    y43m  = (((y400  * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    z43m  = (((z400  * 21.0_DP + c436) * usz + c434) * usz + c432) * uz
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    y43p  =  ((c437  * usy     + c435) * usy + c433) * usy + c431
    z43p  =  ((c437  * usz     + c435) * usz + c433) * usz + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dy43m =  ((dy400 * 21.0_DP + d436) * usy + d434) * usy + c432
    dz43m =  ((dz400 * 21.0_DP + d436) * usz + d434) * usz + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux
    dy43p =  ((d437  * usy     + d435) * usy + d433) * uy
    dz43p =  ((d437  * usz     + d435) * usz + d433) * uz

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    y44m  = (((-y400  * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    z44m  = (((-z400  * 35.0_DP + c446) * usz + c444) * usz + c442) * uz
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    y44p  =  ((c447   * usy     + c445) * usy + c443) * usy + c441
    z44p  =  ((c447   * usz     + c445) * usz + c443) * usz + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dy44m =  ((-dy400 * 35.0_DP + d446) * usy + d444) * usy + c442
    dz44m =  ((-dz400 * 35.0_DP + d446) * usz + d444) * usz + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux
    dy44p =  ((d447   * usy     + d445) * usy + d443) * uy
    dz44p =  ((d447   * usz     + d445) * usz + d443) * uz

    cx1   =  x41p +  x41m
    cx2   =  x42p +  x42m
    cx3   =  x43p +  x43m
    cx4   =  x44p +  x44m
    cx5   =  x44p -  x44m
    cx6   =  x43p -  x43m
    cx7   =  x42p -  x42m
    cx8   =  x41p -  x41m
    cy1   =  y41p +  y41m
    cy2   =  y42p +  y42m
    cy3   =  y43p +  y43m
    cy4   =  y44p +  y44m
    cy5   =  y44p -  y44m
    cy6   =  y43p -  y43m
    cy7   =  y42p -  y42m
    cy8   =  y41p -  y41m
    cz(1) =  z41p +  z41m
    cz(2) =  z42p +  z42m
    cz(3) =  z43p +  z43m
    cz(4) =  z44p +  z44m
    cz(5) =  z44p -  z44m
    cz(6) =  z43p -  z43m
    cz(7) =  z42p -  z42m
    cz(8) =  z41p -  z41m
    dx1   =  dx41p + dx41m
    dx2   =  dx42p + dx42m
    dx3   =  dx43p + dx43m
    dx4   =  dx44p + dx44m
    dx5   =  dx44p - dx44m
    dx6   =  dx43p - dx43m
    dx7   =  dx42p - dx42m
    dx8   =  dx41p - dx41m
    dy1   =  dy41p + dy41m
    dy2   =  dy42p + dy42m
    dy3   =  dy43p + dy43m
    dy4   =  dy44p + dy44m
    dy5   =  dy44p - dy44m
    dy6   =  dy43p - dy43m
    dy7   =  dy42p - dy42m
    dy8   =  dy41p - dy41m
    dz(1) =  dz41p + dz41m
    dz(2) =  dz42p + dz42m
    dz(3) =  dz43p + dz43m
    dz(4) =  dz44p + dz44m
    dz(5) =  dz44p - dz44m
    dz(6) =  dz43p - dz43m
    dz(7) =  dz42p - dz42m
    dz(8) =  dz41p - dz41m


    loop010 : DO l=1,l3d
      w0l =  0.0_DP
      wxl =  0.0_DP
      wyl =  0.0_DP
      wzl =  0.0_DP

      loop020 : DO mz=1,m2
        lz =  iz + mz - 3
        w01 =     cx1 * f3d(l,ix1,iy1,lz) + cx2 * f3d(l,ix2,iy1,lz)       &
          & +     cx3 * f3d(l,ix, iy1,lz) + cx4 * f3d(l,ix4,iy1,lz)       &
          & +     cx5 * f3d(l,ix5,iy1,lz) + cx6 * f3d(l,ix6,iy1,lz)       &
          & +     cx7 * f3d(l,ix7,iy1,lz) + cx8 * f3d(l,ix8,iy1,lz)
        w02 =     cx1 * f3d(l,ix1,iy2,lz) + cx2 * f3d(l,ix2,iy2,lz)       &
          & +     cx3 * f3d(l,ix, iy2,lz) + cx4 * f3d(l,ix4,iy2,lz)       &
          & +     cx5 * f3d(l,ix5,iy2,lz) + cx6 * f3d(l,ix6,iy2,lz)       &
          & +     cx7 * f3d(l,ix7,iy2,lz) + cx8 * f3d(l,ix8,iy2,lz)
        w03 =     cx1 * f3d(l,ix1,iy, lz) + cx2 * f3d(l,ix2,iy, lz)       &
          & +     cx3 * f3d(l,ix, iy, lz) + cx4 * f3d(l,ix4,iy, lz)       &
          & +     cx5 * f3d(l,ix5,iy, lz) + cx6 * f3d(l,ix6,iy, lz)       &
          & +     cx7 * f3d(l,ix7,iy, lz) + cx8 * f3d(l,ix8,iy, lz)
        w04 =     cx1 * f3d(l,ix1,iy4,lz) + cx2 * f3d(l,ix2,iy4,lz)       &
          & +     cx3 * f3d(l,ix, iy4,lz) + cx4 * f3d(l,ix4,iy4,lz)       &
          & +     cx5 * f3d(l,ix5,iy4,lz) + cx6 * f3d(l,ix6,iy4,lz)       &
          & +     cx7 * f3d(l,ix7,iy4,lz) + cx8 * f3d(l,ix8,iy4,lz)
        w05 =     cx1 * f3d(l,ix1,iy5,lz) + cx2 * f3d(l,ix2,iy5,lz)       &
          & +     cx3 * f3d(l,ix, iy5,lz) + cx4 * f3d(l,ix4,iy5,lz)       &
          & +     cx5 * f3d(l,ix5,iy5,lz) + cx6 * f3d(l,ix6,iy5,lz)       &
          & +     cx7 * f3d(l,ix7,iy5,lz) + cx8 * f3d(l,ix8,iy5,lz)
        w06 =     cx1 * f3d(l,ix1,iy6,lz) + cx2 * f3d(l,ix2,iy6,lz)       &
          & +     cx3 * f3d(l,ix, iy6,lz) + cx4 * f3d(l,ix4,iy6,lz)       &
          & +     cx5 * f3d(l,ix5,iy6,lz) + cx6 * f3d(l,ix6,iy6,lz)       &
          & +     cx7 * f3d(l,ix7,iy6,lz) + cx8 * f3d(l,ix8,iy6,lz)
        w07 =     cx1 * f3d(l,ix1,iy7,lz) + cx2 * f3d(l,ix2,iy7,lz)       &
          & +     cx3 * f3d(l,ix, iy7,lz) + cx4 * f3d(l,ix4,iy7,lz)       &
          & +     cx5 * f3d(l,ix5,iy7,lz) + cx6 * f3d(l,ix6,iy7,lz)       &
          & +     cx7 * f3d(l,ix7,iy7,lz) + cx8 * f3d(l,ix8,iy7,lz)
        w08 =     cx1 * f3d(l,ix1,iy8,lz) + cx2 * f3d(l,ix2,iy8,lz)       &
          & +     cx3 * f3d(l,ix, iy8,lz) + cx4 * f3d(l,ix4,iy8,lz)       &
          & +     cx5 * f3d(l,ix5,iy8,lz) + cx6 * f3d(l,ix6,iy8,lz)       &
          & +     cx7 * f3d(l,ix7,iy8,lz) + cx8 * f3d(l,ix8,iy8,lz)

        wxl = ((dx1 * f3d(l,ix1,iy1,lz) + dx2 * f3d(l,ix2,iy1,lz)         &
          & +   dx3 * f3d(l,ix, iy1,lz) + dx4 * f3d(l,ix4,iy1,lz)         &
          & +   dx5 * f3d(l,ix5,iy1,lz) + dx6 * f3d(l,ix6,iy1,lz)         &
          & +   dx7 * f3d(l,ix7,iy1,lz) + dx8 * f3d(l,ix8,iy1,lz)) * cy1  &
          & +  (dx1 * f3d(l,ix1,iy2,lz) + dx2 * f3d(l,ix2,iy2,lz)         &
          & +   dx3 * f3d(l,ix, iy2,lz) + dx4 * f3d(l,ix4,iy2,lz)         &
          & +   dx5 * f3d(l,ix5,iy2,lz) + dx6 * f3d(l,ix6,iy2,lz)         &
          & +   dx7 * f3d(l,ix7,iy2,lz) + dx8 * f3d(l,ix8,iy2,lz)) * cy2  &
          & +  (dx1 * f3d(l,ix1,iy, lz) + dx2 * f3d(l,ix2,iy, lz)         &
          & +   dx3 * f3d(l,ix, iy, lz) + dx4 * f3d(l,ix4,iy, lz)         &
          & +   dx5 * f3d(l,ix5,iy, lz) + dx6 * f3d(l,ix6,iy, lz)         &
          & +   dx7 * f3d(l,ix7,iy, lz) + dx8 * f3d(l,ix8,iy, lz)) * cy3  &
          & +  (dx1 * f3d(l,ix1,iy4,lz) + dx2 * f3d(l,ix2,iy4,lz)         &
          & +   dx3 * f3d(l,ix, iy4,lz) + dx4 * f3d(l,ix4,iy4,lz)         &
          & +   dx5 * f3d(l,ix5,iy4,lz) + dx6 * f3d(l,ix6,iy4,lz)         &
          & +   dx7 * f3d(l,ix7,iy4,lz) + dx8 * f3d(l,ix8,iy4,lz)) * cy4  &
          & +  (dx1 * f3d(l,ix1,iy5,lz) + dx2 * f3d(l,ix2,iy5,lz)         &
          & +   dx3 * f3d(l,ix, iy5,lz) + dx4 * f3d(l,ix4,iy5,lz)         &
          & +   dx5 * f3d(l,ix5,iy5,lz) + dx6 * f3d(l,ix6,iy5,lz)         &
          & +   dx7 * f3d(l,ix7,iy5,lz) + dx8 * f3d(l,ix8,iy5,lz)) * cy5  &
          & +  (dx1 * f3d(l,ix1,iy6,lz) + dx2 * f3d(l,ix2,iy6,lz)         &
          & +   dx3 * f3d(l,ix, iy6,lz) + dx4 * f3d(l,ix4,iy6,lz)         &
          & +   dx5 * f3d(l,ix5,iy6,lz) + dx6 * f3d(l,ix6,iy6,lz)         &
          & +   dx7 * f3d(l,ix7,iy6,lz) + dx8 * f3d(l,ix8,iy6,lz)) * cy6  &
          & +  (dx1 * f3d(l,ix1,iy7,lz) + dx2 * f3d(l,ix2,iy7,lz)         &
          & +   dx3 * f3d(l,ix, iy7,lz) + dx4 * f3d(l,ix4,iy7,lz)         &
          & +   dx5 * f3d(l,ix5,iy7,lz) + dx6 * f3d(l,ix6,iy7,lz)         &
          & +   dx7 * f3d(l,ix7,iy7,lz) + dx8 * f3d(l,ix8,iy7,lz)) * cy7  &
          & +  (dx1 * f3d(l,ix1,iy8,lz) + dx2 * f3d(l,ix2,iy8,lz)         &
          & +   dx3 * f3d(l,ix, iy8,lz) + dx4 * f3d(l,ix4,iy8,lz)         &
          & +   dx5 * f3d(l,ix5,iy8,lz) + dx6 * f3d(l,ix6,iy8,lz)         &
          & +   dx7 * f3d(l,ix7,iy8,lz) + dx8 * f3d(l,ix8,iy8,lz)) * cy8) &
          & *   cz(mz) + wxl

        w00 =  w01 * cy1 + w02 * cy2 + w03 * cy3 + w04 * cy4              &
          & +  w05 * cy5 + w06 * cy6 + w07 * cy7 + w08 * cy8
        w0l =  w00 * cz(mz) + w0l
        wzl =  w00 * dz(mz) + wzl
        wyl = (w01 * dy1 + w02 * dy2 + w03 * dy3 + w04 * dy4              &
          & +  w05 * dy5 + w06 * dy6 + w07 * dy7 + w08 * dy8) * cz(mz)    &
          & +  wyl

      END DO loop020
         
      w0(l) =  w0l
      wx(l) =  wxl / h3x
      wy(l) =  wyl / h3y
      wz(l) =  wzl / h3z

    END DO loop010


    RETURN
  END SUBROUTINE spl3dd

END MODULE spline_mod
