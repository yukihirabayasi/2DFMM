! ************************************************************************** 
     subroutine indus_s(xi,yi,ax,ay)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! Calculate induced velocity at a point from a SOURCE PANEL.                 
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

   ax = 0.0d0
   ay = 0.0d0

   do i = 1, npanel


! ************ normal ******************

    rx = poi(1,1,i) - poi(1,2,i)
    ry = poi(2,1,i) - poi(2,2,i)
    r  = dsqrt(rx*rx + ry*ry)

    gzi_x =  rx/r
    gzi_y =  ry/r
    eta_x = -gzi_y
    eta_y =  gzi_x
    r_j   =  1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)
    r_a   =  0.5d0*ds(i)

    x1 = xi - poi_c(1,i)
    y1 = yi - poi_c(2,i)
    xg = gzi_x*x1 + gzi_y*y1
    yg = eta_x*x1 + eta_y*y1

    if (dabs(yg) .le. 1.0d-20) then
     if (yg .ge. 0.0d0) yg =  1.0d-20
     if (yg .lt. 0.0d0) yg = -1.0d-20
    end if

    coef_a = dlog ((  r_a - xg)*( r_a - xg) + yg*yg) &
   &        -dlog (( -r_a - xg)*(-r_a - xg) + yg*yg)
    coef_b = datan((  r_a - xg) / yg) &
   &        -datan(( -r_a - xg) / yg)

    ax1 = -q(i)*pifour*coef_a
    ay1 =  q(i)*pitwo* coef_b

   if ( dsqrt( x1*x1 + y1*y1 ) .le. 1.0d-10) then
    ax = ax + 0.50d0*q(i)*poi_n(1,i)
    ay = ay + 0.50d0*q(i)*poi_n(2,i)
   else
    ax = ax + ( eta_y*ax1 - gzi_y*ay1)*r_j
    ay = ay + (-eta_x*ax1 + gzi_x*ay1)*r_j
   end if

! ************ mirror (y-direction) ******************

   if (i_mirror .eq. 1) then


    rx = poi(1,1,i) - poi(1,2,i)
    ry = ( 2.0d0*d_mirror - poi(2,1,i) ) - ( 2.0d0*d_mirror - poi(2,2,i) )
    r  = dsqrt(rx*rx + ry*ry)

    gzi_x =  rx/r
    gzi_y =  ry/r
    eta_x = -gzi_y
    eta_y =  gzi_x
    r_j   =  1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)
    r_a   =  0.5d0*ds(i)

    x1 = xi - poi_c(1,i)
    y1 = yi - ( 2.0d0*d_mirror - poi_c(2,i) )
    xg = gzi_x*x1 + gzi_y*y1
    yg = eta_x*x1 + eta_y*y1

    if (dabs(yg) .le. 1.0d-20) then
     if (yg .ge. 0.0d0) yg =  1.0d-20
     if (yg .lt. 0.0d0) yg = -1.0d-20
    end if

    coef_a = dlog ((  r_a - xg)*( r_a - xg) + yg*yg) &
   &        -dlog (( -r_a - xg)*(-r_a - xg) + yg*yg)
    coef_b = datan((  r_a - xg) / yg) &
   &        -datan(( -r_a - xg) / yg)

    ax1 = -q(i)*pifour*coef_a
    ay1 =  q(i)*pitwo* coef_b
    ax = ax + ( eta_y*ax1 - gzi_y*ay1)*r_j
    ay = ay + (-eta_x*ax1 + gzi_x*ay1)*r_j


   end if

! ************ end mirror image ******************


   end do

   return
   end

