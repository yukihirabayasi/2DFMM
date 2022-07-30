! ************************************************************************** 
     subroutine biot_b(xi,yi,ax,ay)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Calculate induced velocity at a point(x,y) from a vortex blob element     
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   ax = 0.0d0
   ay = 0.0d0

   do i = 1, nvor_b

! ************** normal **************

     rx = xi - vor_b(1,i)
     ry = yi - vor_b(2,i)
     r = dsqrt( rx**2 + ry**2 )
     rv2 = 1.0d0 / ( r**2 )

     if ( r .gt. 1.0d-6 ) then
       xai = r / cor_b(i)
       rv = cir_b(i) * ( 1.0d0 - dexp( -xai**2 ) )
       ax = ax - ry*rv*rv2
       ay = ay + rx*rv*rv2
     end if

! ************** mirror (y-direction) **************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     rx = xi - vor_b(1,i)
     ry = yi - ( 2.0d0*d_mirror - vor_b(2,i) )
     r = dsqrt( rx**2 + ry**2 )

     if ( r .gt. 1.0d-6 ) then
       rv2 =  1.0d0 / r**2
       xai =  r / cor_b(i)
       rv  = -cir_b(i) * ( 1.0d0 - dexp( -xai**2 ) )
       ax  =  ax - ry*rv*rv2
       ay  =  ay + rx*rv*rv2
     end if

! ==========
     end if 
! ==========

! ************** end mirror image **************

   end do

   ax = ax * pitwo
   ay = ay * pitwo

  return
  end


! ************************************************************************** 
     subroutine biot_s(xi,yi,ax,ay)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Calculate induced velocity at a point(x,y) from a vortex sheet element    
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   ax = 0.0d0
   ay = 0.0d0

   do i = 1, nvor_s

! ************** normal **************

     rx = vor_s(1,1,i) - vor_s(1,2,i)
     ry = vor_s(2,1,i) - vor_s(2,2,i)
     r = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x
     r_j   =  1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

     r_a = 0.5d0 * r
     r_b = cor_s(i)
     r_s = ( 2.0d0*r_a )*( 2.0d0*r_b )
     omg = cir_s(i) / r_s

     x1 = xi - vor_sc(1,i)
     y1 = yi - vor_sc(2,i)
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

     ax = ax + (  eta_y*ax1 - gzi_y*ay1 )*r_j
     ay = ay + ( -eta_x*ax1 + gzi_x*ay1 )*r_j

! ************ mirror (y-direction) ******************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     rx = vor_s(1,1,i) - vor_s(1,2,i)
     ry = ( 2.0d0*d_mirror - vor_s(2,1,i) ) - ( 2.0d0*d_mirror - vor_s(2,2,i) )
     r  = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x = rx*r_inv
     gzi_y = ry*r_inv
     eta_x =-gzi_y
     eta_y = gzi_x
     r_j   = 1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )

     r_a =  0.5d0 * r
     r_b =  cor_s(i)
     r_s =  ( 2.0d0*r_a )*( 2.0d0*r_b )
     omg = -cir_s(i) / r_s

     x1 = xi - vor_sc(1,i)
     y1 = yi - ( 2.0d0*d_mirror - vor_sc(2,i) )
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

     ax = ax + (  eta_y*ax1 - gzi_y*ay1 ) * r_j
     ay = ay + ( -eta_x*ax1 + gzi_x*ay1 ) * r_j

! ==========
     end if 
! ==========

! ************ end mirror image ******************

   end do

  return
  end


