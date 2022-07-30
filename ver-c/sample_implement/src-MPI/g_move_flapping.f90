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

   time_in_cycle = time - cycle_time * int(time/cycle_time)

   write(*,*) ''
   write(*,*) '  time              =', time
   write(*,*) '  time step         =', nt
   write(*,*) '  time / cycle_time =', time/cycle_time


!--------------------
   do i = 1 , npanel 
!--------------------

   k = panel_id(i)

!--- heaving location ---

    dx = 0.0d0

   if ( ( time_in_cycle .ge. 0.0d0 ) .and. &
        ( time_in_cycle .le. t_t ) ) then

    dy = - time_in_cycle

   else if ( ( time_in_cycle .gt. t_t )         .and. &
             ( time_in_cycle .le. t_t + dt_t ) ) then

    dy = - t_t - dt_t/pi * &
                 (  dsin( pi * ( time_in_cycle - 2.0d0*t_t ) / dt_t) &
                  + dsin( pi * t_t / dt_t ) )

   else if ( ( time_in_cycle .gt. t_t + dt_t )               .and. &
             ( time_in_cycle .le. 0.5d0*cycle_time + t_t ) )  then

    dy = time_in_cycle - 2.0d0*t_t - dt_t &
                       - dt_t/pi * (  dsin ( pi*( dt_t - t_t ) / dt_t ) &
                                    + dsin ( pi*t_t / dt_t ) )

   else if ( ( time_in_cycle .gt. 0.5d0*cycle_time + t_t )          .and. &
             ( time_in_cycle .le. 0.5d0*cycle_time + t_t + dt_t ) )  then

    dy = 0.5d0*cycle_time - t_t - dt_t &
        - dt_t/pi &
      * (  dsin( pi*( dt_t - t_t) / dt_t) &
         - dsin( pi*( time_in_cycle - 0.5d0*cycle_time - 2.0d0*t_t ) / dt_t ) )

   else if ( ( time_in_cycle .gt. 0.5d0*cycle_time + t_t + dt_t ) .and. &
             ( time_in_cycle .le. 1.0d0*cycle_time ) )             then

    dy = cycle_time - time_in_cycle

   end if


!--- pitching angle ---

   if ( ( time_in_cycle .ge. 0.0d0 ) .and. &
        ( time_in_cycle .le. t_r ) ) then

    theta_p = pitch_angle_base_down

   else if ( ( time_in_cycle .gt. t_r )         .and. &
             ( time_in_cycle .le. t_r + dt_r ) ) then

    theta_p =  pitch_angle_base_down &
             + 0.50d0 * d_alpha_0 * ( time_in_cycle - t_r ) &
             - 0.25d0 * d_alpha_0 / pi * dt_r &
               * ( dsin( 2.0d0 * pi * (time_in_cycle - 2.0d0 * t_r ) / dt_r ) &
                  +dsin( 2.0d0 * pi * t_r / dt_r) )

   else if ( ( time_in_cycle .gt. t_r + dt_r )               .and. &
             ( time_in_cycle .le. 0.5d0*cycle_time + t_r ) )  then

    theta_p = pitch_angle_base_up

   else if ( ( time_in_cycle .gt. 0.5d0*cycle_time + t_r )          .and. &
             ( time_in_cycle .le. 0.5d0*cycle_time + t_r + dt_r ) )  then

    theta_p =  pitch_angle_base_up &
             - 0.5d0*d_alpha_0*( time_in_cycle - 0.5d0*cycle_time - 1.0d0*t_r)&
             + 0.25d0 * d_alpha_0 / pi * dt_r &
               * (  dsin ( 2.0d0 * pi * (  time_in_cycle &
                                         - 0.5d0*cycle_time -2.0d0*t_r)/dt_r) &
                  + dsin ( 2.0d0 * pi * t_r / dt_r ) )

   else if ( ( time_in_cycle .gt. 0.5d0*cycle_time + t_r + dt_r) .and. &
             ( time_in_cycle .le. 1.0d0*cycle_time ) )            then

    theta_p = pitch_angle_base_down

   end if


    parts_position(1,k) = dx
    parts_position(2,k) = dy
    parts_position(3,k) = theta_p*180.0d0/pi

!---------------
   do j = 1 , 2 
!---------------

     poi_old (1,j,i) = poi  (1,j,i)
     poi_old (2,j,i) = poi  (2,j,i)

     poi     (1,j,i) =   poi_b (1,j,i)*dcos( theta_p ) &
                       + poi_b (2,j,i)*dsin( theta_p )
     poi     (2,j,i) = - poi_b (1,j,i)*dsin( theta_p ) &
                       + poi_b (2,j,i)*dcos( theta_p )

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


