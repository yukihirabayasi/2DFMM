! ************************************************************************** 
     subroutine read_file
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

   open (11,file='../cal_cond/cond.dat',status='old')
    read(11,*) theta
    read(11,*) omz
    read(11,*) i_continue
    read(11,*) re
    read(11,*) dt
    read(11,*) cut_r
    read(11,*) i_panel_height_check
    read(11,*) c_h
    read(11,*) loop
    read(11,*) i_write
    read(11,*) file1
    read(11,*) file2
    read(11,*) i_change_cond
    read(11,*) i_kutta_cond
    read(11,*) i_mirror
    read(11,*) d_mirror
    read(11,*) i_edge
    read(11,*) i_layer
    read(11,*) 
    read(11,*) 
    read(11,*) r_rot
    read(11,*) tsr
    read(11,*) 
    read(11,*) 
    read(11,*) cycle_time
    read(11,*) pitch_angle_base_up
    read(11,*) pitch_angle_base_down
    read(11,*) dt_t
    read(11,*) dt_r
    read(11,*) dt_lag
    read(11,*) 
    read(11,*) 
    read(11,*) velocity
    read(11,*) accelerate
    read(11,*) uniform_velocity
    read(11,*) u_accelerate
    read(11,*) 
    read(11,*) 
    read(11,*) freq
    read(11,*) amp
    read(11,*) 
    read(11,*) 
    read(11,*) k_p
    read(11,*) s_int
    read(11,*) s_amp
    read(11,*) f_cx
    read(11,*) f_cy
   close(11)

! ////////// Uniform Velocity //////////*

   theta0 = theta * pi / 180.0d0

   uinf = uniform_velocity
   vinf = 0.0d0

   velocity = -velocity
   velocity_base = velocity

   uniform_velocity_base = uniform_velocity
   dt_base = dt

! ////////// for wind turbine ///////////

   v_rot = tsr*dsqrt( uinf**2 + vinf**2 )

! ////////// for flapping ///////////

   pitch_angle_base_up   = pitch_angle_base_up   * pi / 180.0d0
   pitch_angle_base_down = pitch_angle_base_down * pi / 180.0d0
   dt_t   = dt_t   * cycle_time
   dt_r   = dt_r   * cycle_time
   dt_lag = dt_lag * cycle_time
   t_t = 0.25d0*cycle_time - 0.5d0*dt_t
   t_r = 0.25d0*cycle_time - 0.5d0*dt_r - 1.0d0*dt_lag
   d_alpha_0 = 2.0d0*( pitch_angle_base_up - pitch_angle_base_down ) / dt_r

   write(*,*)                  '  '
   write(*,*)                  '  '
   write(*,*)                    '*=========================================*'
   write(*,*)                    '*  Numerical Simulation of Unsteady Flow  *'
   write(*,*)                    '*  around bodies by a 2-D Vortex Method   *'
   write(*,*)                    '*    for Parallel Processor with MPI      *'
   write(*,*)                    '*    coded by Kota Fukuda (2004/03/14)    *'
   write(*,*)                    '*=========================================*'
   write(*,*)                  '  '
   write(*,*)                  '  '
   write(*,*)                  '*** Calculation condition ***'
   write(*,'(1x, A, f10.2)')   ' Reynolds number   =', re
   write(*,'(1x, A, f10.5)')   ' Time interval     =', dt
   write(*,'(1x, A, f10.5)')   ' Attack angle      =', theta
   write(*,'(1x, A, f10.5)')   ' Calculation space =', cut_r
   write(*,'(1x, A, f10.5)')   ' Uniform flow (x)  =', uinf
   write(*,'(1x, A, f10.5)')   ' Uniform flow (y)  =', vinf
   write(*,'(1x, A, f10.5)')   ' Angular velocity  =', omz

  return
  end


! ************************************************************************** 
     subroutine read_change_file
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

   open (11,file='../cal_cond/cond.dat',status='old')
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) re
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) r_rot
    read(11,*) tsr
    read(11,*) 
    read(11,*) 
    read(11,*) cycle_time
    read(11,*) pitch_angle_base_up
    read(11,*) pitch_angle_base_down
    read(11,*) dt_t
    read(11,*) dt_r
    read(11,*) dt_lag
    read(11,*) 
    read(11,*) 
    read(11,*) 
    read(11,*) accelerate
    read(11,*) 
    read(11,*) u_accelerate
   close(11)

   write(*,*)                  '  '
   write(*,*)                  '  '
   write(*,*)                    '*=========================================*'
   write(*,*)                    '*  Numerical Simulation of Unsteady Flow  *'
   write(*,*)                    '*  around bodies by a 2-D Vortex Method   *'
   write(*,*)                    '*    for Parallel Processor with MPI      *'
   write(*,*)                    '*    coded by Kota Fukuda (2004/03/14)    *'
   write(*,*)                    '*=========================================*'
   write(*,*)                  '  '
   write(*,*)                  '  '
   write(*,*)                  '*** Calculation condition ***'
   write(*,'(1x, A, f10.2)')   ' Reynolds number   =', re
   write(*,'(1x, A, f10.5)')   ' Time interval     =', dt
   write(*,'(1x, A, f10.5)')   ' Attack angle      =', theta
   write(*,'(1x, A, f10.5)')   ' Calculation space =', cut_r
   write(*,'(1x, A, f10.5)')   ' Uniform flow (x)  =', uinf
   write(*,'(1x, A, f10.5)')   ' Uniform flow (y)  =', vinf
   write(*,'(1x, A, f10.5)')   ' Angular velocity  =', omz

    return
    end

