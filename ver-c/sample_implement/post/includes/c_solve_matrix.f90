! ************************************************************************** 
     subroutine solve_matrix                                                 
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

   allocate( a(npanel, npanel),b(npanel),ipiv(npanel) )

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

     rx   = poi(1,1,j) - poi(1,2,j)
     ry   = poi(2,1,j) - poi(2,2,j)
     r    = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x

     r_j  = 1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)
     r_a  = 0.5d0 * ds(j)

     x1   = poi_c(1,i) - poi_c(1,j)
     y1   = poi_c(2,i) - poi_c(2,j)
     xg   = gzi_x*x1 + gzi_y*y1
     yg   = eta_x*x1 + eta_y*y1

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

     if ( dsqrt( x1**2 + y1**2 ) .le. 1.0d-10 ) then
       ax = 0.5d0*poi_n(1,i)
       ay = 0.5d0*poi_n(2,i)
     end if

! ************** mirror (y-direction) **************

! ===============================
     if ( i_mirror .eq. 1 ) then 
! ===============================

     rx   = poi(1,1,j) - poi(1,2,j)
     ry   = ( 2.0d0*d_mirror - poi(2,1,j) ) - ( 2.0d0*d_mirror - poi(2,2,j) )
     r    = dsqrt( rx**2 + ry**2 )
     r_inv = 1.0d0 / r

     gzi_x =  rx * r_inv
     gzi_y =  ry * r_inv
     eta_x = -gzi_y
     eta_y =  gzi_x

     r_j  = 1.0d0 / (gzi_x*eta_y - gzi_y*eta_x)
     r_a  = 0.5d0 * ds(j)

     x1   = poi_c(1,i) - poi_c(1,j)
     y1   = poi_c(2,i) - ( 2.0d0*d_mirror - poi_c(2,j) )
     xg   = gzi_x*x1 + gzi_y*y1
     yg   = eta_x*x1 + eta_y*y1

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

!     if ( dsqrt( x1**2 + y1**2 ) .le. 1.0d-10 ) then
!       ax = ax + 0.5d0*poi_n(1,i)
!       ay = ay - 0.5d0*poi_n(2,i)
!     end if

! ==========
     end if 
! ==========

! ************** end mirror image **************

     a(i,j) = ax*poi_n(1,i) + ay*poi_n(2,i)

   end do
   end do

! ************************************ 
! ---// set boundary conditionin //--- 
! ************************************ 

   b(1:npanel) = 0.0d0

   do i = 1, npanel
     x1 = poi_c(1,i)
     y1 = poi_c(2,i)
     call indus_v( x1,y1,u2,v2 )
     call biot_b ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )
     vmx  = uinf + (-omz*y1) - vm(1,1,i) + ( u2 + u3 + u4 )
     vmy  = vinf + ( omz*x1) - vm(2,1,i) + ( v2 + v3 + v4 )
     b(i) = -( vmx*poi_n(1,i) + vmy*poi_n(2,i) )
   end do

! ************************* 
! ---// solve matrix  //--- 
! ************************* 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ---------------// DGESV //-----------------------
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   call dgesv( npanel,1,a,npanel,ipiv,b,npanel,info )

   do i = 1, npanel
    q(i) = b(i)
   end do

!   write(*,'(A29,i8)') ' lapack dgesv info =         ',info

   if ( info .ne. 0 ) then
     write(*,*) 'Velocity matrix solver error !!'
     stop
   end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   deallocate( a,b,ipiv )

  return
  end


! ************************************************************************** 
     subroutine non_slip_condition_check
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!     Check the velocity on the boundary point : uw(i,j)                     
!      uw(a,b)*n(a,b) .eq. 0.0   ==> ok                                      
!                     .ne. 0.0   ==> no                                      
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   ncount = 0

   do i = 1, npanel

     x1 = poi_c(1,i)
     y1 = poi_c(2,i)
     call indus_s( x1,y1,u1,v1 )
     call indus_v( x1,y1,u2,v2 )
     call biot_b ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )
     vmx = uinf + (-omz*y1) - vm(1,1,i)
     vmy = vinf + ( omz*x1) - vm(2,1,i)
     ux  = vmx + ( u1 + u2 + u3 + u4 )
     uy  = vmy + ( v1 + v2 + v3 + v4 )

     v_nor = ux*poi_n(1,i) + uy*poi_n(2,i)
     diff  = dabs( v_nor )

     if ( diff .gt. 1.0d-10 ) then
       ncount = ncount + 1
       write(*,*) ' error! ( not b.c. ) '
       write(*,*) ' i = ',i
       write(*,*) ' differnt of vn=', v_nor
       stop
     end if

   end do

   if ( ncount .ne. 0 ) write(*,*) 'error count = ',ncount,'/',npanel


!///////////// Cal. Velocity on Vortex Panel ////////////

   do i = 1, npanel
     x1 = poi_c(1,i) + dh * poi_n(1,i)
     y1 = poi_c(2,i) + dh * poi_n(2,i)
     call indus_s( x1,y1,u1,v1 )
     call indus_v( x1,y1,u2,v2 )
     call biot_b ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )
     vmx = uinf + (-omz*y1) - vm(1,1,i)
     vmy = vinf + ( omz*x1) - vm(2,1,i)
     uw(1,i) = vmx + ( u1 + u2 + u3 + u4 )
     uw(2,i) = vmy + ( v1 + v2 + v3 + v4 )
   end do

  return
  end


