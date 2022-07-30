! ************************************************************************** 
     subroutine solve_matrix_pressure                                        
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

   external dgesv

   integer,          allocatable :: ipiv(:)
   double precision, allocatable :: a (:,:)
   double precision, allocatable :: b   (:)

   allocate( a(npanel,npanel),b(npanel),ipiv(npanel) )

! ************************** 
! ---// initialization //--- 
! ************************** 

   a(1:npanel,1:npanel) = 0.0d0

! ********************** 
! ---// set matrix //--- 
! ********************** 

   do i = 1, npanel
   do j = 1, npanel

     ax   = 0.0d0
     ay   = 0.0d0

! ************** normal **************

     rx = poi(1,1,j) - poi(1,2,j)
     ry = poi(2,1,j) - poi(2,2,j)
     r  = dsqrt( rx**2 + ry**2 )
     r_inv =  1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x

     r_j   = 1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )
     r_a   = 0.5d0 * ds(j)

! ////////// Derection of vector (x1,y1) is important !! //////////*

     x1 = poi_c(1,j) - poi_c(1,i)
     y1 = poi_c(2,j) - poi_c(2,i)
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     if ( dabs(yg) .le. 1.0d-20 ) then
       if ( yg .ge. 0.0d0 ) yg =  1.0d-20
       if ( yg .lt. 0.0d0 ) yg = -1.0d-20
     end if

     rxg1 = xg + r_a
     rxg2 = xg - r_a

     coef_a = dlog( rxg1**2 + yg**2 ) - dlog( rxg2**2 + yg**2 )
     coef_b = datan( rxg1 / yg ) - datan( rxg2 / yg )

     ax1 = pifour*coef_a
     ay1 = pitwo *coef_b

     ax  = ax + (  eta_y*ax1 - gzi_y*ay1 )*r_j
     ay  = ay + ( -eta_x*ax1 + gzi_x*ay1 )*r_j

! ************** mirror (y-direction) **************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     rx = poi(1,1,j) - poi(1,2,j)
     ry = ( 2.0d0*d_mirror - poi(2,1,j) ) - ( 2.0d0*d_mirror - poi(2,2,j) )
     r  = dsqrt( rx**2 + ry**2 )
     r_inv =  1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x

     r_j   = 1.0d0 / ( gzi_x*eta_y - gzi_y*eta_x )
     r_a   = 0.5d0 * ds(j)

! ////////// Derection of vector (x1,y1) is important !! //////////*

     x1 = poi_c(1,j) - poi_c(1,i)
     y1 = ( 2.0d0*d_mirror - poi_c(2,j) ) - poi_c(2,i)
     xg = gzi_x*x1 + gzi_y*y1
     yg = eta_x*x1 + eta_y*y1

     if ( dabs(yg) .le. 1.0d-20 ) then
       if ( yg .ge. 0.0d0 ) yg =  1.0d-20
       if ( yg .lt. 0.0d0 ) yg = -1.0d-20
     end if

     rxg1 = xg + r_a
     rxg2 = xg - r_a

     coef_a = dlog( rxg1**2 + yg**2 ) - dlog( rxg2**2 + yg**2 )
     coef_b = datan( rxg1 / yg ) - datan( rxg2 / yg )

     ax1 = pifour*coef_a
     ay1 = pitwo *coef_b

     ax  = ax + (  eta_y*ax1 - gzi_y*ay1 )*r_j
     ay  = ay + ( -eta_x*ax1 + gzi_x*ay1 )*r_j

! ==========
     end if 
! ==========

! ************ end mirror image ******************

     a(i,j) = ax*poi_n(1,j) + ay*poi_n(2,j)

   end do
   end do

! ************************************ 
! ---// set boundary conditionin //--- 
! ************************************ 

   icf = min(1, nt)

   b(1:npanel) = 0.0d0

   do i = 1, npanel
     x1 = poi_c(1,i)
     y1 = poi_c(2,i)
     call gh_b( x1,y1,i,bh_b )
     call gh_s( x1,y1,i,bh_s )
     call gh_r( x1,y1,i,bh_r )
     call gh_a( x1,y1,i,bh_a )
     b(i) = -( bh_b + bh_s + bh_a + icf*bh_r )
   end do

! ************************* 
! ---// solve matrix  //--- 
! ************************* 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ---------------// DGESV //-----------------------
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   call dgesv( npanel,1,a,npanel,ipiv,b,npanel,info )

   do i = 1, npanel
     ph(i) = b(i)
   end do

!   write(*,'(A29,i8)') ' lapack dgesv info =         ',info

   if ( info .ne. 0 ) then
     write(*,*) 'Pressure matrix solver error !!'
     stop
   end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   deallocate( a,b,ipiv )

  return
  end


