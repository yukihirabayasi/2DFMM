! ************************************************************************** 
     subroutine drift_b
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   new_b = 0

   do i = 1, nvor_b

     x0 = vor_b(1,i)
     y0 = vor_b(2,i)

     cir2 = cir_b(i)
     cor2 = cor_b(i)

!/////////// modified by K.Kuji /2015/10/15/ ////////////
    
!//////////// calculate by 3rd Adams-Bashforth ///////////

     vmx = 23.0d0/12.0d0*vor_vb(1,1,i) &
    &    -  4.0d0/ 3.0d0*vor_vb(1,2,i) &
    &    +  5.0d0/12.0d0*vor_vb(1,3,i) 

     vmy = 23.0d0/12.0d0*vor_vb(2,1,i) &
    &    -  4.0d0/ 3.0d0*vor_vb(2,2,i) &
    &    +  5.0d0/12.0d0*vor_vb(2,3,i)

!//////////// calculate by 4th Runge-Kutta ///////////

! <<<<< Notyet! Cording !! >>>>>

!/////////////////// end modified ///////////////////////

     x1 = x0 + vmx*dt
     y1 = y0 + vmy*dt

     if ( i_mirror .eq. 1 .and. y1 .le.  d_mirror ) y1 = 2.0d0*d_mirror - y1

     call cross2( x1,y1,cir2,cor2,ic,height,icross )

     r1 = dsqrt( x1**2 + y1**2 )

     if ( icross .eq. 1 ) cycle

! =================================================
     if ( icross .eq. 0 .and. r1 .le. cut_r ) then 
! =================================================
     new_b = new_b + 1
     vor_b   (1,new_b) = x1
     vor_b   (2,new_b) = y1
     cir_b     (new_b) = cir2
     cor_b     (new_b) = cor2
     vor_vb(1,1,new_b) = vor_vb(1,1,i)
     vor_vb(2,1,new_b) = vor_vb(2,1,i)
     vor_vb(1,2,new_b) = vor_vb(1,2,i)
     vor_vb(2,2,new_b) = vor_vb(2,2,i)
     vor_vb(1,3,new_b) = vor_vb(1,3,i)
     vor_vb(2,3,new_b) = vor_vb(2,3,i)
     blob_id   (new_b) = blob_id   (i)
! ========
     else 
! ========
     do j = 1, nw
       if ( blob_id(i) .eq. j ) cir_r(j) = cir_r(j) + cir2
     end do
! ==========
     end if 
! ==========

   end do

   nvor_b = new_b

  return
  end


! ************************************************************************** 
      subroutine drift_s
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   new_b = 0
   new_s = 0

   r_unit = 2.00d0 * dh

   do i = 1, nvor_s

     x0 = vor_sc(1,i)
     y0 = vor_sc(2,i)

     cor2 = cor_s(i)
     cir2 = cir_s(i)

!/////////// modified by K.Kuji /2015/10/15/ //////////// 

!//////////// calculate by 3rd Adams-Bashforth /////////// 

     vmx = 23.0d0/12.0d0*vor_vs(1,1,i) &
    &    -  4.0d0/ 3.0d0*vor_vs(1,2,i) &
    &    +  5.0d0/12.0d0*vor_vs(1,3,i) 

     vmy = 23.0d0/12.0d0*vor_vs(2,1,i) &
    &    -  4.0d0/ 3.0d0*vor_vs(2,2,i) &
    &    +  5.0d0/12.0d0*vor_vs(2,3,i)

!//////////// calculate by 4th Runge-Kutta ///////////

! <<<<< Notyet! Cording !! >>>>>

!/////////////////// end modified ///////////////////////

     x1 = x0 + vmx*dt
     y1 = y0 + vmy*dt

     if ( i_mirror .eq. 1 .and. y1 .le.  d_mirror ) y1 = 2.0d0*d_mirror - y1

     call cross2( x1,y1,cir2,cor2,ic,height,icross )

     if ( icross .eq. 1 ) cycle

     if ( icross .eq. 0 ) then

       c1 = 2.0d0 * cor2
       c2 = vor_sr(i)

! ====================================================
       if ( height .ge. r_unit .or. c1 .ge. c2 ) then 
! ====================================================

       n_new = idint( c2 / c1 )
       if ( n_new .eq. 0 ) n_new = 1

       c2 = vor_sr(i) / dble( n_new )

       do k = 1, n_new
         new_b = new_b + 1
         vor_nb   (1,new_b) = ( x1-0.5d0*dble(n_new-1)*c2 ) + dble(k-1)*c2
         vor_nb   (2,new_b) = y1
         vor_nvb(1,1,new_b) = vor_vs(1,1,i)
         vor_nvb(2,1,new_b) = vor_vs(2,1,i)
         vor_nvb(1,2,new_b) = vor_vs(1,2,i)
         vor_nvb(2,2,new_b) = vor_vs(2,2,i)
         vor_nvb(1,3,new_b) = vor_vs(1,3,i)
         vor_nvb(2,3,new_b) = vor_vs(2,3,i)
         cir_nb     (new_b) = cir2 / dble( n_new )
         cor_nb     (new_b) = dsqrt( 2.0d0 * cor2 * c2 / pi )
         blob_nid   (new_b) = sheet_id  (i)
       end do

! ==========
       else 
! ==========

       new_s = new_s + 1
       gx = poi(1,1,ic) - poi(1,2,ic)
       gy = poi(2,1,ic) - poi(2,2,ic)
       gr = dsqrt( gx**2 + gy**2 )
       vor_sc  (1,new_s) = x1
       vor_sc  (2,new_s) = y1
       vor_s (1,1,new_s) = vor_sc(1,new_s) + 0.5d0 * vor_sr(i) * gx/gr
       vor_s (2,1,new_s) = vor_sc(2,new_s) + 0.5d0 * vor_sr(i) * gy/gr
       vor_s (1,2,new_s) = vor_sc(1,new_s) - 0.5d0 * vor_sr(i) * gx/gr
       vor_s (2,2,new_s) = vor_sc(2,new_s) - 0.5d0 * vor_sr(i) * gy/gr
       vor_vs(1,1,new_s) = vor_vs(1,1,i)
       vor_vs(2,1,new_s) = vor_vs(2,1,i)
       vor_vs(1,2,new_s) = vor_vs(1,2,i)
       vor_vs(2,2,new_s) = vor_vs(2,2,i)
       vor_vs(1,3,new_s) = vor_vs(1,3,i)
       vor_vs(2,3,new_s) = vor_vs(2,3,i)
       cir_s     (new_s) = cir2
       cor_s     (new_s) = cor2
       vor_sr    (new_s) = vor_sr    (i)
       sheet_id  (new_s) = sheet_id  (i)

! ============
       end if 
! ============
     end if

   end do

   nvor_s = new_s

   do i = 1, new_b
     nvor_b = nvor_b + 1
     vor_b   (1,nvor_b) = vor_nb   (1,i)
     vor_b   (2,nvor_b) = vor_nb   (2,i)
     vor_vb(1,1,nvor_b) = vor_nvb(1,1,i)
     vor_vb(2,1,nvor_b) = vor_nvb(2,1,i)
     vor_vb(1,2,nvor_b) = vor_nvb(1,2,i)
     vor_vb(2,2,nvor_b) = vor_nvb(2,2,i)
     vor_vb(1,3,nvor_b) = vor_nvb(1,3,i)
     vor_vb(2,3,nvor_b) = vor_nvb(2,3,i)
     cir_b     (nvor_b) = cir_nb     (i)
     cor_b     (nvor_b) = cor_nb     (i)
     blob_id   (nvor_b) = blob_nid   (i)
   end do

  return
  end


! ************************************************************************** 
      subroutine cross2( xpc,ypc,cir2,cor2,ic,height,icross )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! DISPOSITION OF VORTEX ELEMENTS ON WALL VICINITY ( REFLECTIONÅEUPTHRUST )   
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   icross = 0

   dmin = 1.0d10

   do i = 1, npanel

    poicx = poi_c(1,i) + dh*poi_n(1,i)
    poicy = poi_c(2,i) + dh*poi_n(2,i)
    dmin2 = dsqrt( (xpc - poicx)**2 + (ypc - poicy)**2 )

    if ( dmin2 .lt. dmin ) then
      dmin = dmin2
      ic = i
    end if

   end do

   ax1 = ( poi(1,1,ic) + dh*poi_n(1,ic) ) - xpc
   ay1 = ( poi(2,1,ic) + dh*poi_n(2,ic) ) - ypc
   ax2 = ( poi(1,2,ic) + dh*poi_n(1,ic) ) - xpc
   ay2 = ( poi(2,2,ic) + dh*poi_n(2,ic) ) - ypc
   ss  = 0.5d0 * ( ax1*ay2 - ay1*ax2 )
   hei = dabs( ss * 2.0d0 / ds(ic) )

   poinx  = poi_n(1,ic)
   poiny  = poi_n(2,ic)
   poicx  = poi_c(1,ic) + dh*poi_n(1,ic)
   poicy  = poi_c(2,ic) + dh*poi_n(2,ic)
   vss    = dsqrt( (xpc - poicx)**2 + (ypc - poicy)**2 )
   vin    = poinx*( xpc - poicx ) + poiny*( ypc - poicy )
   height = hei

   rx1 = poi(1,1,ic) - poicx
   ry1 = poi(2,1,ic) - poicy
   r1  = dsqrt( rx1**2 + ry1**2 )

   rx2 = poi(1,2,ic) - poicx
   ry2 = poi(2,2,ic) - poicy
   r2  = dsqrt( rx2**2 + ry2**2 )

   fg1 = ( ax1*rx1 + ay1*ry1 ) / r1
   fg2 = ( ax2*rx2 + ay2*ry2 ) / r2

   if ( fg1 .lt. 0.0d0 .or. fg2 .lt. 0.0d0 ) goto 100

! ====================================
   if ( ( vin/vss ) .ge. 0.0d0 ) then 
! ====================================
   if ( hei .lt. cor2 ) then
     cir2 = 0.5d0 * ( cor2 + hei ) * cir2 / cor2
     cor2 = 0.5d0 * ( cor2 + hei )
     xpc = xpc + poinx * ( cor2 - hei )
     ypc = ypc + poiny * ( cor2 - hei )
     if ( dabs(cir2) .lt. 1.0d-8 ) icross = 1
     if ( cor2       .lt. 1.0d-8 ) icross = 1
   end if
! ======
   else 
! ======
   if ( hei .lt. cor2 ) then
     cir2 = 0.5d0 * ( cor2 - hei ) * cir2 / cor2
     cor2 = 0.5d0 * ( cor2 - hei )
     xpc  = xpc + poinx * ( cor2 + hei )
     ypc  = ypc + poiny * ( cor2 + hei )
     if ( dabs(cir2) .lt. 1.0d-8 ) icross = 1
     if ( cor2       .lt. 1.0d-8 ) icross = 1
   else
     icross = 1
   end if
! ========
   end if 
! ========

   100 continue

  return
  end


