! ************************************************************************** 
     subroutine gh_b(xi,yi,ip,bh)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Calculate induced pressure at a point(x,y) from a vortex BLOB element     
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   bh = 0.0d0

   do i = 1, nvor_b

! ************ normal ******************

! ////////// Derection of vector R is important !! //////////*

     rx = vor_b(1,i) - xi
     ry = vor_b(2,i) - yi
     r  = dsqrt( rx**2 + ry**2 )
     rv2 = 1.0d0 / ( r**2 )

     if( r .gt. 1.0d-6 ) then
       xai =  r / cor_b(i)
       rv  =  cir_b(i) * ( 1.0d0 - dexp( -xai**2 ) )
       ax  = -ry*rv*rv2
       ay  =  rx*rv*rv2
       bh  = bh + pitwo*( ax*( vor_vb(1,1,i) - vm(1,1,ip) ) &
      &                 + ay*( vor_vb(2,1,i) - vm(2,1,ip) ) )
     end if

! ************ mirror (y-direction) ************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

! ////////// Derection of vector R is important !! //////////*

     rx = vor_b(1,i) - xi
     ry = ( 2.0d0*d_mirror - vor_b(2,i) ) - yi
     rv2 = 1.0d0 / ( r**2 )

     if ( r .gt. 1.0d-6 ) then
       xai =  r/cor_b(i)
       rv  = -cir_b(i) * ( 1.0d0 - dexp( -xai**2 ) )
       ax  = -ry*rv*rv2
       ay  =  rx*rv*rv2
       bh  = bh + pitwo*( ax*(  vor_vb(1,1,i) - vm(1,1,ip) ) &
      &                 + ay*( -vor_vb(2,1,i) + vm(2,1,ip) ) )
     end if

! ==========
     end if 
! ==========

! ************ end mirror image ************

   end do

  return
  end


! ************************************************************************** 
     subroutine gh_s(xi,yi,ip,bh)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Calculate induced pressure at a point(x,y) from a vortex SHEET element    
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   bh = 0.0d0

   do i = 1, nvor_s

! ************** normal **************

     rx = vor_s(1,1,i) - vor_s(1,2,i)
     ry = vor_s(2,1,i) - vor_s(2,2,i)
     r  = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r 

     gzi_x =  rx*r_inv
     gzi_y =  ry*r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x
     r_j   =  1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

     r_a = 0.5d0 * r
     r_b = cor_s(i)
     r_s = ( 2.0d0*r_a )*( 2.0d0*r_b )
     omg = cir_s(i) / r_s

! ////////// Derection of vector (x1,y1) is important !! //////////*

     x1 = vor_sc(1,i) - xi
     y1 = vor_sc(2,i) - yi
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     r_ax = xg + r_a
     r_bx = xg - r_a
     r_cy = yg + r_b
     r_dy = yg - r_b

     if ( dabs(r_ax) .le. 1.0d-20 ) then
       if ( r_ax .ge. 0.0d0 ) r_ax =  1.0d-20
       if ( r_ax .lt. 0.0d0 ) r_ax = -1.0d-20
     end if

     if ( dabs(r_bx) .le. 1.0d-20 ) then
       if ( r_bx .ge. 0.0d0 ) r_bx =  1.0d-20
       if ( r_bx .lt. 0.0d0 ) r_bx = -1.0d-20
     end if

     if ( dabs(r_cy) .le. 1.0d-20 ) then
       if ( r_cy .ge. 0.0d0 ) r_cy =  1.0d-20
       if ( r_cy .lt. 0.0d0 ) r_cy = -1.0d-20
     end if

     if ( dabs(r_dy) .le. 1.0d-20 ) then
       if ( r_dy .ge. 0.0d0 ) r_dy =  1.0d-20
       if ( r_dy .lt. 0.0d0 ) r_dy = -1.0d-20
     end if

     ce_a1 =  r_dy*datan( r_ax/r_dy ) + 0.5d0*r_ax*dlog( r_ax**2 + r_dy**2 )
     ce_a2 = -r_cy*datan( r_ax/r_cy ) - 0.5d0*r_ax*dlog( r_ax**2 + r_cy**2 )
     ce_a3 = -r_dy*datan( r_bx/r_dy ) - 0.5d0*r_bx*dlog( r_bx**2 + r_dy**2 )
     ce_a4 =  r_cy*datan( r_bx/r_cy ) + 0.5d0*r_bx*dlog( r_bx**2 + r_cy**2 )

     ce_b1 =  r_bx*datan( r_dy/r_bx ) + 0.5d0*r_dy*dlog( r_dy**2 + r_bx**2 )
     ce_b2 = -r_ax*datan( r_dy/r_ax ) - 0.5d0*r_dy*dlog( r_dy**2 + r_ax**2 )
     ce_b3 = -r_bx*datan( r_cy/r_bx ) - 0.5d0*r_cy*dlog( r_cy**2 + r_bx**2 )
     ce_b4 =  r_ax*datan( r_cy/r_ax ) + 0.5d0*r_cy*dlog( r_cy**2 + r_ax**2 )

     ax1 = omg * pitwo * ( ce_a1 + ce_a2 + ce_a3 + ce_a4 )
     ay1 = omg * pitwo * ( ce_b1 + ce_b2 + ce_b3 + ce_b4 )

     ax = (  eta_y*ax1 - gzi_y*ay1 )*r_j
     ay = ( -eta_x*ax1 + gzi_x*ay1 )*r_j

     bh = bh + ( ax*( vor_vs(1,1,i) - vm(1,1,ip) ) &
    &          + ay*( vor_vs(2,1,i) - vm(2,1,ip) ) )

! ************** mirror (y-direction) **************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     rx = vor_s(1,1,i) - vor_s(1,2,i)
     ry = ( 2.0d0*d_mirror - vor_s(2,1,i) ) - ( 2.0d0*d_mirror - vor_s(2,2,i) )
     r = dsqrt(rx*rx + ry*ry)
     r_inv = 1.0d0 / r 

     gzi_x =  rx*r_inv
     gzi_y =  ry*r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x
     r_j   =  1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

     r_a =  0.5d0 * r
     r_b =  cor_s(i)
     r_s =  ( 2.0d0*r_a )*( 2.0d0*r_b )
     omg = -cir_s(i) / r_s

! ////////// Derection of vector R is important !! //////////*

     x1 = vor_sc(1,i) - xi
     y1 = ( 2.0d0*d_mirror - vor_sc(2,i) ) - yi
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     r_ax = xg + r_a
     r_bx = xg - r_a
     r_cy = yg + r_b
     r_dy = yg - r_b

     if ( dabs(r_ax) .le. 1.0d-20 ) then
       if ( r_ax .ge. 0.0d0 ) r_ax =  1.0d-20
       if ( r_ax .lt. 0.0d0 ) r_ax = -1.0d-20
     end if

     if ( dabs(r_bx) .le. 1.0d-20 ) then
       if ( r_bx .ge. 0.0d0 ) r_bx =  1.0d-20
       if ( r_bx .lt. 0.0d0 ) r_bx = -1.0d-20
     end if

     if ( dabs(r_cy) .le. 1.0d-20 ) then
       if ( r_cy .ge. 0.0d0 ) r_cy =  1.0d-20
       if ( r_cy .lt. 0.0d0 ) r_cy = -1.0d-20
     end if

     if ( dabs(r_dy) .le. 1.0d-20 ) then
       if ( r_dy .ge. 0.0d0 ) r_dy =  1.0d-20
       if ( r_dy .lt. 0.0d0 ) r_dy = -1.0d-20
     end if

     ce_a1 =  r_dy*datan( r_ax/r_dy ) + 0.5d0*r_ax*dlog( r_ax**2 + r_dy**2 )
     ce_a2 = -r_cy*datan( r_ax/r_cy ) - 0.5d0*r_ax*dlog( r_ax**2 + r_cy**2 )
     ce_a3 = -r_dy*datan( r_bx/r_dy ) - 0.5d0*r_bx*dlog( r_bx**2 + r_dy**2 )
     ce_a4 =  r_cy*datan( r_bx/r_cy ) + 0.5d0*r_bx*dlog( r_bx**2 + r_cy**2 )

     ce_b1 =  r_bx*datan( r_dy/r_bx ) + 0.5d0*r_dy*dlog( r_dy**2 + r_bx**2 )
     ce_b2 = -r_ax*datan( r_dy/r_ax ) - 0.5d0*r_dy*dlog( r_dy**2 + r_ax**2 )
     ce_b3 = -r_bx*datan( r_cy/r_bx ) - 0.5d0*r_cy*dlog( r_cy**2 + r_bx**2 )
     ce_b4 =  r_ax*datan( r_cy/r_ax ) + 0.5d0*r_cy*dlog( r_cy**2 + r_ax**2 )

     ax1 = omg * pitwo * ( ce_a1 + ce_a2 + ce_a3 + ce_a4 )
     ay1 = omg * pitwo * ( ce_b1 + ce_b2 + ce_b3 + ce_b4 )

     ax = (  eta_y*ax1 - gzi_y*ay1 )*r_j
     ay = ( -eta_x*ax1 + gzi_x*ay1 )*r_j

     bh = bh + ( ax*(  vor_vs(1,1,i) - vm(1,1,ip) ) &
    &          + ay*( -vor_vs(2,1,i) + vm(2,1,ip) ) )

! ==========
     end if 
! ==========

! ************** end mirror image **************

   end do

  return
  end


! ************************************************************************** 
     subroutine gh_r(xi,yi,ip,bh)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! INDUCED PRESSURE BY VISCOUS EFFECT                                         
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   bh = 0.0d0

   do i = 1, npanel

! ************** normal **************

     ! j = 1
       ux = lay_vc(1,1,i)
       uy = lay_vc(2,1,i)
       vtt = ( uy*poi_n(1,i) - ux*poi_n(2,i) )

     if ( i_layer .eq. 1 ) then
       do j = 1, nlayer(i)
         ux = lay_vc(1,j,i) - lay_vc(1,j-1,i)
         uy = lay_vc(2,j,i) - lay_vc(2,j-1,i)
         vtt = vtt + ( uy*poi_n(1,i) - ux*poi_n(2,i) )
       end do
     end if

     omega_s = vtt / dh

     b1 = omega_s
     b2 = omega_s

     rx = poi(1,1,i) - poi(1,2,i)
     ry = poi(2,1,i) - poi(2,2,i)
     r  = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x
     r_j   =  1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

! ////////// Derection of vector (x1,y1) is important !! //////////*

     x1 = xi - poi_c(1,i)
     y1 = yi - poi_c(2,i)
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     if ( dabs(yg) .le. 1.0d-20 ) then
       if ( yg .ge. 0.0d0 ) yg =  1.0d-20
       if ( yg .lt. 0.0d0 ) yg = -1.0d-20
     end if

     r_a     = 0.5d0 * ds(i)
     r_a_1_4 = 1.0d0 / ( 4.0d0 * r_a )

     rxg1 =  xg + r_a
     rxg2 =  xg - r_a

     coef_a = dlog( rxg1**2 + yg**2 ) - dlog( rxg2**2 + yg**2 )
     coef_b = datan( rxg1 / yg ) - datan( rxg2 / yg )

     ax1 = ( -coef_a*yg + 2.0d0*coef_b*xg + 2.0d0*r_a*coef_b )       * r_a_1_4
     ax2 = (  coef_a*yg - 2.0d0*coef_b*xg + 2.0d0*r_a*coef_b )       * r_a_1_4
     ay1 = ( -coef_a*xg - 2.0d0*coef_b*yg + 4.0d0*r_a - r_a*coef_a ) * r_a_1_4
     ay2 = (  coef_a*xg + 2.0d0*coef_b*yg - 4.0d0*r_a - r_a*coef_a ) * r_a_1_4

     axx1 = (  eta_y*ax1 - gzi_y*ay1 )*r_j
     axx2 = (  eta_y*ax2 - gzi_y*ay2 )*r_j
     ayy1 = ( -eta_x*ax1 + gzi_x*ay1 )*r_j
     ayy2 = ( -eta_x*ax2 + gzi_x*ay2 )*r_j

     axx1 = axx1 * pitwo * b1
     axx2 = axx2 * pitwo * b2
     ayy1 = ayy1 * pitwo * b1
     ayy2 = ayy2 * pitwo * b2

     ax  = axx1 + axx2
     ay  = ayy1 + ayy2

     bh = bh + ( poi_n(1,i)*ax + poi_n(2,i)*ay ) / re

! ************** mirror (y-direction) **************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     ! j = 1
       ux = lay_vc(1,1,i)
       uy = lay_vc(2,1,i)
       vtt = ( uy*poi_n(1,i) - ux*poi_n(2,i) )

     if ( i_layer .eq. 1 ) then
       do j = 1, nlayer(i)
         ux = lay_vc(1,j,i) - lay_vc(1,j-1,i)
         uy = lay_vc(2,j,i) - lay_vc(2,j-1,i)
         vtt = vtt + ( uy*poi_n(1,i) - ux*poi_n(2,i) )
       end do
     end if

     omega_s = vtt / dh

     b1 = -omega_s
     b2 = -omega_s

     rx = poi(1,1,i) - poi(1,2,i)
     ry = ( 2.0d0*d_mirror - poi(2,1,i) ) - ( 2.0d0*d_mirror - poi(2,2,i) )
     r  = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x
     r_j   =  1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

! ////////// Derection of vector (x1,y1) is important !! //////////*

     x1 = xi - poi_c(1,i)
     y1 = yi - ( 2.0d0*d_mirror - poi_c(2,i) )
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     if ( dabs(yg) .le. 1.0d-20 ) then
       if ( yg .ge. 0.0d0 ) yg =  1.0d-20
       if ( yg .lt. 0.0d0 ) yg = -1.0d-20
     end if

     r_a     = 0.5d0 * ds(i)
     r_a_1_4 = 1.0d0 / ( 4.0d0*r_a )

     rxg1 =  xg + r_a
     rxg2 =  xg - r_a

     coef_a = dlog( rxg1**2 + yg**2 ) - dlog( rxg2**2 + yg**2 )
     coef_b = datan( rxg1 / yg ) - datan( rxg2 / yg )

     ax1 = ( -coef_a*yg + 2.0d0*coef_b*xg + 2.0d0*r_a*coef_b )       * r_a_1_4
     ax2 = (  coef_a*yg - 2.0d0*coef_b*xg + 2.0d0*r_a*coef_b )       * r_a_1_4
     ay1 = ( -coef_a*xg - 2.0d0*coef_b*yg + 4.0d0*r_a - r_a*coef_a ) * r_a_1_4
     ay2 = (  coef_a*xg + 2.0d0*coef_b*yg - 4.0d0*r_a - r_a*coef_a ) * r_a_1_4

     axx1 = (  eta_y*ax1 - gzi_y*ay1 )*r_j
     axx2 = (  eta_y*ax2 - gzi_y*ay2 )*r_j
     ayy1 = ( -eta_x*ax1 + gzi_x*ay1 )*r_j
     ayy2 = ( -eta_x*ax2 + gzi_x*ay2 )*r_j

     axx1 = axx1 * pitwo * b1
     axx2 = axx2 * pitwo * b2
     ayy1 = ayy1 * pitwo * b1
     ayy2 = ayy2 * pitwo * b2

     ax  = axx1 + axx2
     ay  = ayy1 + ayy2

     bh = bh + ( poi_n(1,i)*ax + poi_n(2,i)*ay ) / re

! ==========
     end if 
! ==========

! ************** end mirror image **************

   end do

  return
  end


! ************************************************************************** 
     subroutine gh_a(xi,yi,ip,bh)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! INDUCED PRESSURE BY INERTIAL FORCE EFFECT                                  
!============================================================================

  use global
  use parameter
  use math_constant
  use gauss_legendre
  implicit none
  include 'inc_2d'

   bh = 0.0d0

   call intp_data

   do i = 1, npanel

! ************** normal **************

   do j = 1, n_weight

     xx = 0.5d0*( 1.0d0 - phai(j) )*poi(1,1,i) &
    &   + 0.5d0*( 1.0d0 + phai(j) )*poi(1,2,i)
     yy = 0.5d0*( 1.0d0 - phai(j) )*poi(2,1,i) &
    &   + 0.5d0*( 1.0d0 + phai(j) )*poi(2,2,i)

! ////////// Derection of vector R is important !! //////////*

     rx = xx - xi
     ry = yy - yi
     rc = dsqrt( rx**2 + ry**2 )

     if ( rc .lt. 1.0d-6 ) cycle

     ax = poi_n(1,i)*am(1,1,i) + poi_n(2,i)*am(2,1,i)
     bh = bh + ( pitwo * ax * dlog(rc) * ds(i) * 0.5d0 * w(j) )

   end do

! ************** mirror (y-direction) **************

! =============================
   if ( i_mirror .eq. 1 ) then 
! =============================

   do j = 1, n_weight

     xx = 0.5d0*( 1.0d0 - phai(j) )*poi(1,1,i) &
    &   + 0.5d0*( 1.0d0 + phai(j) )*poi(1,2,i)
     yy = 0.5d0*( 1.0d0 - phai(j) )*( 2.0d0*d_mirror - poi(2,1,i) ) &
    &   + 0.5d0*( 1.0d0 + phai(j) )*( 2.0d0*d_mirror - poi(2,2,i) )

! ////////// Derection of vector R is important !! //////////*

     rx = xx - xi
     ry = yy - yi
     rc = dsqrt( rx**2 + ry**2 )

     if ( rc .lt. 1.0d-6 ) cycle

     ax = poi_n(1,i)*am(1,1,i) + poi_n(2,i)*am(2,1,i)
     bh = bh + ( pitwo * ax * dlog(rc) * ds(i) * 0.5d0 * w(j) )

   end do

! ========
   end if 
! ========

! ************** end mirror image **************

   end do

  return
  end


