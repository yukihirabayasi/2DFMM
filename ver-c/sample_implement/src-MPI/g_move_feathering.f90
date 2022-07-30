! ************************************************************************** 
      subroutine move
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

  double precision, allocatable :: poi_old(:,:,:)

   allocate ( poi_old (2,2,npanel) )

! ***** FETHERING ANGLE *****

   if ( time .eq. 0.0d0 ) then
     s_old = 0.0d0
   end if

   omega_f = 2.0d0 * pi * k_p * time
   s       = ( s_int + s_amp * dsin( omega_f ) ) * pi / 180.0d0
   ds_p    = s - s_old
   s_old   = s

   write(*,*)             ' '
   write(*,'(A18,f12.5)') 'feathering angle =', s*180.0d0/pi

! ***************************

!---------------
   do k = 1, nw 
!---------------

    if ( time .eq. 0.0d0 ) then
      parts_position(1,k) = 0.0d0
      parts_position(2,k) = 0.0d0
    end if

     parts_position(1,k) = parts_position(1,k)
     parts_position(2,k) = parts_position(2,k)
     parts_position(3,k) = 0.0d0

!--------
  end do 
!--------

!--------------------
   do i = 1 , npanel 
!--------------------

!---------------
   do j = 1 , 2 
!---------------

    if ( time .eq. 0.0 ) then
      poi (1,j,i) = poi_b (1,j,i)
      poi (2,j,i) = poi_b (2,j,i)
    end if

    poi_old (1,j,i) = poi  (1,j,i)
    poi_old (2,j,i) = poi  (2,j,i)
    poi     (1,j,i) = ( poi_old(1,j,i) - f_cx ) * dcos(ds_p) + ( poi_old(2,j,i) - f_cy ) * dsin( ds_p ) + f_cx
    poi     (2,j,i) = ( poi_old(2,j,i) - f_cy ) * dcos(ds_p) - ( poi_old(1,j,i) - f_cx ) * dsin( ds_p ) + f_cy

!--------------
   end do
!--------------

    panel_c_x     = 0.5d0 * ( poi     (1,1,i) + poi     (1,2,i) )
    panel_c_y     = 0.5d0 * ( poi     (2,1,i) + poi     (2,2,i) )
    panel_c_x_old = 0.5d0 * ( poi_old (1,1,i) + poi_old (1,2,i) )
    panel_c_y_old = 0.5d0 * ( poi_old (2,1,i) + poi_old (2,2,i) )

    vm (1,2,i) = vm (1,1,i)
    vm (2,2,i) = vm (2,1,i)
    am (1,2,i) = am (1,1,i)
    am (2,2,i) = am (2,1,i)

    vm (1,1,i) = ( panel_c_x - panel_c_x_old ) / dt
    vm (2,1,i) = ( panel_c_y - panel_c_y_old ) / dt
    am (1,1,i) = ( vm(1,1,i) - vm(1,2,i) )     / dt
    am (2,1,i) = ( vm(2,1,i) - vm(2,2,i) )     / dt

    if ( time .eq. 0.0d0 ) then
      vm (1,1,i) = 0.0d0
      vm (2,1,i) = 0.0d0
      am (1,1,i) = 0.0d0
      am (2,1,i) = 0.0d0
    end if

!--------------
   end do
!--------------

   deallocate( poi_old )

   return

   end

