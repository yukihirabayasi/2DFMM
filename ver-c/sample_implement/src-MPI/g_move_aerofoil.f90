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

   allocate( poi_old (2,2,npanel) )

!   velocity = velocity - 0.5d0 * accelerate * dt
   velocity = velocity - accelerate * dt	 !modify by kuji /2014/12/10

   dx = velocity * dt
   dy = 0.0d0

!   uniform_velocity = uniform_velocity - 0.5d0 * u_accelerate * dt
   uniform_velocity = uniform_velocity + u_accelerate * dt 	! modify by kuji /2014/12/10

   uinf = uniform_velocity   ! add by K.F. (2012/08/13)
   vinf = 0.0d0              ! add by K.F. (2012/08/13)

!===============
   do k = 1, nw 
!===============

   if ( time .eq. 0.0 ) then
     parts_position(1,k) = 0.0d0
     parts_position(2,k) = 0.0d0
   end if

   parts_position(1,k) = parts_position(1,k) + dx
   parts_position(2,k) = parts_position(2,k) + dy
   parts_position(3,k) = 0.0d0

!=========
   end do 
!=========

   do i = 1 , npanel 

     do j = 1 , 2 

       if ( time .eq. 0.0 ) then
         poi (1,j,i) = poi_b (1,j,i)
         poi (2,j,i) = poi_b (2,j,i)
       end if

       poi_old (1,j,i) = poi (1,j,i)
       poi_old (2,j,i) = poi (2,j,i)
       poi     (1,j,i) = poi (1,j,i) + dx
       poi     (2,j,i) = poi (2,j,i) + dy

     end do

     panel_c_x     = 0.5d0 * ( poi(1,1,i) + poi(1,2,i) )
     panel_c_y     = 0.5d0 * ( poi(2,1,i) + poi(2,2,i) )
     panel_c_x_old = 0.5d0 * ( poi_old(1,1,i) + poi_old(1,2,i) )
     panel_c_y_old = 0.5d0 * ( poi_old(2,1,i) + poi_old(2,2,i) )

     vm (1,2,i) = vm (1,1,i)
     vm (2,2,i) = vm (2,1,i)
     am (1,2,i) = am (1,1,i)
     am (2,2,i) = am (2,1,i)

     vm (1,1,i) = ( panel_c_x - panel_c_x_old ) / dt
     vm (2,1,i) = ( panel_c_y - panel_c_y_old ) / dt
     am (1,1,i) = ( vm(1,1,i) - vm(1,2,i) ) / dt
     am (2,1,i) = ( vm(2,1,i) - vm(2,2,i) ) / dt

     if ( time .eq. 0.0 ) then
       vm (1,1,i) = 0.0d0
       vm (2,1,i) = 0.0d0
       am (1,1,i) = 0.0d0
       am (2,1,i) = 0.0d0
     end if

   end do

   deallocate( poi_old )

  return
  end


