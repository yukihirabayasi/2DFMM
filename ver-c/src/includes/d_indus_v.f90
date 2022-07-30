! ************************************************************************** 
   subroutine indus_v(xi,yi,ax,ay)
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! Calculate induced velocity at a point from a VORTEX PANEL.                 
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

    b1 = cir_g(i)
    b2 = cir_g(i)

    rx = poi(1,1,i) - poi(1,2,i)
    ry = poi(2,1,i) - poi(2,2,i)
    r  = dsqrt(rx*rx + ry*ry)

    gzi_x =  rx/r
    gzi_y =  ry/r
    eta_x = -gzi_y
    eta_y =  gzi_x
    r_j   =  1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)

    x1 = xi - poi_c(1,i)
    y1 = yi - poi_c(2,i)
    xg = gzi_x*x1 + gzi_y*y1
    yg = eta_x*x1 + eta_y*y1

    if (dabs(yg) .le. 1.0d-20) then
     if (yg .ge. 0.0d0) yg =  1.0d-20
     if (yg .lt. 0.0d0) yg = -1.0d-20
    end if

    r_a = 0.5d0*ds(i)
    r_a_1_4 = 1.0d0 / (4.0d0*r_a)

    coef_a = dlog ( (r_a - xg)**2 + yg**2) &
   &        -dlog ( (r_a + xg)**2 + yg**2)
    coef_b = datan( (r_a - xg)/yg) &
   &        -datan((-r_a - xg)/yg)

    ax1  = ( yg*coef_a + 2.0d0*xg*coef_b + 2.0d0*r_a*coef_b)*r_a_1_4
    ax2  = (-yg*coef_a - 2.0d0*xg*coef_b + 2.0d0*r_a*coef_b)*r_a_1_4
    ay1  = ( 4.0d0*r_a + coef_a*xg + coef_a*r_a - 2.0d0*coef_b*yg)*r_a_1_4
    ay2  = (-4.0d0*r_a - coef_a*xg + coef_a*r_a + 2.0d0*coef_b*yg)*r_a_1_4
    axx1 = ( eta_y*ax1 - gzi_y*ay1)*r_j
    axx2 = ( eta_y*ax2 - gzi_y*ay2)*r_j
    ayy1 = (-eta_x*ax1 + gzi_x*ay1)*r_j
    ayy2 = (-eta_x*ax2 + gzi_x*ay2)*r_j

    axx1 = axx1*pitwo*b1
    axx2 = axx2*pitwo*b2
    ayy1 = ayy1*pitwo*b1
    ayy2 = ayy2*pitwo*b2

    ax = ax + axx1 + axx2
    ay = ay + ayy1 + ayy2


! ************ mirror (y-direction) ******************

   if (i_mirror .eq. 1) then


    b1 = - cir_g(i)
    b2 = - cir_g(i)

    rx = poi(1,1,i) - poi(1,2,i)
    ry = ( 2.0d0*d_mirror - poi(2,1,i) ) - ( 2.0d0*d_mirror - poi(2,2,i) )
    r  = dsqrt(rx*rx + ry*ry)

    gzi_x =  rx/r
    gzi_y =  ry/r
    eta_x = -gzi_y
    eta_y =  gzi_x
    r_j   =  1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)

    x1 = xi - poi_c(1,i)
    y1 = yi - ( 2.0d0*d_mirror - poi_c(2,i) )
    xg = gzi_x*x1 + gzi_y*y1
    yg = eta_x*x1 + eta_y*y1

    if (dabs(yg) .le. 1.0d-20) then
     if (yg .ge. 0.0d0) yg =  1.0d-20
     if (yg .lt. 0.0d0) yg = -1.0d-20
    end if

    r_a = 0.5d0*ds(i)
    r_a_1_4 = 1.0d0 / (4.0d0*r_a)

    coef_a = dlog ( (r_a - xg)**2 + yg**2) &
   &        -dlog ( (r_a + xg)**2 + yg**2)
    coef_b = datan( (r_a - xg)/yg) &
   &        -datan((-r_a - xg)/yg)

    ax1  = ( yg*coef_a + 2.0d0*xg*coef_b + 2.0d0*r_a*coef_b)*r_a_1_4
    ax2  = (-yg*coef_a - 2.0d0*xg*coef_b + 2.0d0*r_a*coef_b)*r_a_1_4
    ay1  = ( 4.0d0*r_a + coef_a*xg + coef_a*r_a - 2.0d0*coef_b*yg)*r_a_1_4
    ay2  = (-4.0d0*r_a - coef_a*xg + coef_a*r_a + 2.0d0*coef_b*yg)*r_a_1_4
    axx1 = ( eta_y*ax1 - gzi_y*ay1)*r_j
    axx2 = ( eta_y*ax2 - gzi_y*ay2)*r_j
    ayy1 = (-eta_x*ax1 + gzi_x*ay1)*r_j
    ayy2 = (-eta_x*ax2 + gzi_x*ay2)*r_j

    axx1 = axx1*pitwo*b1
    axx2 = axx2*pitwo*b2
    ayy1 = ayy1*pitwo*b1
    ayy2 = ayy2*pitwo*b2

    ax = ax + axx1 + axx2
    ay = ay + ayy1 + ayy2


   end if

! ************ end mirror image ******************


   end do

   return
   end


