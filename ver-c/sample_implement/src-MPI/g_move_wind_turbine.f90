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

 double precision, dimension(:,:,:), allocatable :: poi_old

  allocate ( poi_old (2,2,npanel) )

!--------------------
   do i = 1 , npanel 
!--------------------

   k = panel_id(i)

   theta_p = 2.0d0*pi/dble(nw)*dble(k-1) + v_rot/r_rot*time

     dx = r_rot*dcos(theta_p)
     dy = r_rot*dsin(theta_p)

    parts_position(1,k) = dx
    parts_position(2,k) = dy
    parts_position(3,k) = (-theta_p+0.5d0*pi)/pi*180.0d0

!---------------
   do j = 1 , 2 
!---------------

     poi_old (1,j,i) = poi  (1,j,i)
     poi_old (2,j,i) = poi  (2,j,i)

     poi     (1,j,i) =   poi_b (1,j,i)*dcos( - theta_p + 0.5d0*pi ) &
                       + poi_b (2,j,i)*dsin( - theta_p + 0.5d0*pi )
     poi     (2,j,i) = - poi_b (1,j,i)*dsin( - theta_p + 0.5d0*pi ) &
                       + poi_b (2,j,i)*dcos( - theta_p + 0.5d0*pi )

     poi     (1,j,i) = poi (1,j,i) + dx
     poi     (2,j,i) = poi (2,j,i) + dy

!--------------
   end do
!--------------

    panel_c_x     = 0.5d0*( poi     (1,1,i) + poi     (1,2,i) )
    panel_c_y     = 0.5d0*( poi     (2,1,i) + poi     (2,2,i) )

    panel_c_x_old = 0.5d0*( poi_old (1,1,i) + poi_old (1,2,i) )
    panel_c_y_old = 0.5d0*( poi_old (2,1,i) + poi_old (2,2,i) )

    vm (1,2,i) = vm (1,1,i)
    vm (2,2,i) = vm (2,1,i)

    am (1,2,i) = am (1,1,i)
    am (2,2,i) = am (2,1,i)

    vm (1,1,i) = ( panel_c_x - panel_c_x_old ) / dt
    vm (2,1,i) = ( panel_c_y - panel_c_y_old ) / dt

    am (1,1,i) = ( vm(1,1,i) - vm(1,2,i) ) / dt
    am (2,1,i) = ( vm(2,1,i) - vm(2,2,i) ) / dt

    if (time .eq. 0.0) then
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


