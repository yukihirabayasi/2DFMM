! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                    
!     Title:		2-D Vortical Flow Solver                     
!     Purpose:		Numerical Simulation of Unsteady Flow        
!			through a Pump by a 2-D Vortex Method        
!                                                                    
!     Coded by:    	Kota Fukuda    (2004/03/14)                  
!     Modified by: 	Wataru Hattori (2019/08/15)                  
!                                                                    
!     Place:       	Yokohama National University, Japan          
!                       University of Maryland, USA                  
!			JAXA / JEDI center, Japan                    
!			Tokai University, Japan                      
!                                                                    
!     Version: 		ver.003-005                                  
!                                                                    
!     Description:                                                   
!     -----------                                                    
!       Code type:		Original code                        
!       Vortex model:		blob type, sheet type                
!       Viscous diffusion:	Core spreading model                 
!       Time development :	3rd order Adams-Bashforth Method     
!                                                                    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                    
!  Copyright:                                                        
!  ---------                                                         
!                                                                    
!   This program was written by Kota Fukuda, Tokai University, Japan 
!                                                                    
!                                                                    
!   This is UNPUBLISHED PROPRIETARY SOURCE CODE:                     
!   -------------------------------------------                      
!                                                                    
!        The contents of this file may not be disclosed,             
!        copied or duplicated in any form, in whole or in part,      
!        without the prior written permission of Kota Fukuda.        
!                                                                    
!        Redistribution in source and binary forms,                  
!        with or without modification, are prohibited                
!        except redistribution among permitted members.              
!                                                                    
!        The programs are available on "as is" basis.                
!                                                                    
!        They are not guaranteed to be free of bugs,                 
!        and the author assumes no responsibility                    
!        for any potential problems.                                 
!                                                                    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  include 'a_modules.f90'

!*****************************
!***** // Main Program // ****
!*****************************

  use global
  implicit none
  include 'inc_2d'
  include 'mpif.h'

   call system_clock( count_rate = hz )

!  ////////// mpi initialize ///////////

   if ( num_procs .ne. 1 ) call mpi_initial

!  ////////// process check ///////////

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )
   if ( num_procs .ne. 1 ) then
     if ( my_rank .eq. 0 ) write(*,*) '  '
     if ( my_rank .eq. 0 ) write(*,*) '  '
     if ( my_rank .eq. 0 ) write(*,*) '***** Process check *****'
   call MPI_BARRIER( MPI_COMM_WORLD,ierr )
     write(*,'(2x,A,i4,x,A,i4,A)') 'Process',my_rank,' /',num_procs,'  is alive'
   end if

!  ////////// input cal. condition and panel data ///////////

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

!=============================
   if ( my_rank .eq. 0 ) then 
!=============================

   call read_file

   if ( i_continue .eq. 1 ) then
     call continued_file
     if( i_change_cond .eq. 1 ) call read_change_file
   else
     nt     = 0
     nvor_b = 0
     nvor_s = 0
     time   = 0.0d0
     call body
     call move
     call panel
     call panel_height
   end if

   call file_open

!  ////////// check of License expiration date etc ///////////

   call initial_check
   call create_run_check_file

!=========
   end if 
!=========

!  ////////// bcast data from RANK = 0 ///////////

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )
   if ( num_procs .ne. 1 ) then
     call mpi_bcast_cal_data
     call mpi_bcast_panel_data
     call mpi_bcast_vor_data
   end if
   call MPI_BARRIER( MPI_COMM_WORLD,ierr )


! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
!  ---------------------------  
!  ---// START MAIN LOOP //---  
!  ---------------------------  
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

   10 continue

!=============================
   if ( my_rank .eq. 0 ) then 
!=============================

!  ////////// check step //////////

   total_time = MPI_WTIME()

   write(*,*) ''
   write(*,*) ''
   write(*,*) '-------',nt,'step            -------'

!  ////////// data output //////////

   if ( mod(nt,i_write) .eq. 0 ) call data_out

   call run_check_and_close_run

!  ////////// solve velocity matrix ///////////

   call system_clock( count = clock_1 )

  time_matrix_velo = MPI_WTIME()

   call solve_matrix
   call non_slip_condition_check

  time_matrix_velo = MPI_WTIME() - time_matrix_velo

   call run_check_and_close_run

!  ////////// cal. yplus ///////////

   call layer

!  ////////// velo data update ///////////
   
   if ( i_method .eq. 0 ) call velo_update
 
!=========
   end if 
!=========

!  ////////// bcast data from RANK = 0 ///////////

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )
  time_trans = MPI_WTIME()
   if ( num_procs .ne. 1 ) then
     call mpi_bcast_cal_update_data
     call mpi_bcast_panel_update_data
     call mpi_bcast_vor_data
   end if
  time_trans = MPI_WTIME() - time_trans
   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

!  ////////// calculation of velocity field //////////

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )
  time_get_velo = MPI_WTIME()
   if ( num_procs .ne. 1 ) then
     call get_velo_mpi
!     call get_velo
   end if
  time_get_velo = MPI_WTIME() - time_get_velo
   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

!=============================
   if ( my_rank .eq. 0 ) then 
!=============================

   call run_check_and_close_run

!  ////////// solve pressure matrix ///////////

  time_matrix_pres = MPI_WTIME()
   call solve_matrix_pressure
  time_matrix_pres = MPI_WTIME() - time_matrix_pres

   call run_check_and_close_run

!  ////////// calculation of force ///////////

   call force

!  ////////// time integration //////////

   nt   = nt + 1
   time = time+dt

!  ////////// viscous diffusion //////////

  time_core = MPI_WTIME()
   call core
  time_core = MPI_WTIME() - time_core

   call run_check_and_close_run

!  ////////// panel deformation and move ///////////

   call move
   call panel

!  ////////// drift vortex elements ///////////

  time_drift = MPI_WTIME()
   call drift_b
   call drift_s
  time_drift = MPI_WTIME() - time_drift

   call run_check_and_close_run

!  ////////// introduction of vortex elements ///////////

  time_nascent = MPI_WTIME()
   call nascent
   if ( i_kutta_cond .eq. 1 ) call kutta_condition
   call kelvine
  time_nascent = MPI_WTIME() - time_nascent

   call run_check_and_close_run

!  ////////// time step check ///////////

   if ( nt .gt. loop ) goto 40

   fquit = 0
   call MPI_BCAST( fquit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )

!  ////////// check total time ///////////

   total_time = MPI_WTIME() - total_time 

!  ////////// data output //////////

   call file_write

   write(*,*) ''
   write(*,*) '  nvor_b,     nvor_s'
   write(*,'(i9,3x,i9)')  nvor_b, nvor_s
   write(*,*) ''
   write(*,*) '  time step,    time '
   write(*,'(i9,3x,f9.5)') nt, time

!   dt = dt_base * (uniform_velocity_base - velocity_base) &
!                 /(uniform_velocity      - velocity)

   goto 10

! >>>>>>>>>>>>>>>>>>>>>>>>>>> 
! --------------------------- 
!  ---// END MAIN LOOP //---  
! --------------------------- 
! >>>>>>>>>>>>>>>>>>>>>>>>>>> 


!==============================
   else ! ( rank not equal 0 ) 
!==============================

   call MPI_BCAST( fquit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
   if ( fquit .eq. 1 ) then 
     write(*,*) 'quit RANK: ', my_rank
     goto 100
   end if

   goto 10

!=========
   end if 
!=========

  40 continue

   call file_close

  100 continue

   if ( my_rank .eq. 0 .and. num_procs .ne. 1 ) then
     fquit = 1
     call MPI_BCAST( fquit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
   end if

!  ////////// mpi finalize ///////////

   call MPI_FINALIZE( ierr )

  end



! >>>>>>>>>>>>>>>>>>>>> 
! *** include files *** 
! >>>>>>>>>>>>>>>>>>>>> 

! ////////// Subroutine A //////////

   include './a_run_control.f90'

! ////////// Subroutine B //////////

   include './b_read_cond_file.f90'
   include './b_continued_file.f90'
   include './b_initial_check.f90'
   include './b_panel.f90'
   include './b_body.f90'

! ////////// Subroutine C //////////

   include './c_solve_matrix.f90'
   include './c_solve_matrix_pressure.f90'

! ////////// Subroutine D //////////

   include './d_biot.f90'
   include './d_indus.f90'
   include './d_biot_p.f90'
!   include './d_tree.f90'

! ////////// Subroutine E //////////

   include './e_core.f90'
   include './e_drift.f90'
   include './e_nascent.f90'
   include './e_velo.f90'
   include './e_velo_mpi.f90'
   include './e_velo_update.f90'
   include './e_layer.f90'

! ////////// Subroutine F //////////

   include './f_data_out.f90'
   include './f_file.f90'
   include './f_force.f90'

! ////////// Subroutine G //////////

   include './g_kutta.f90'
   include './g_move_aerofoil.f90'
!   include './g_move_feathering.f90'
!   include './g_non-move.f90'
!   include './g_move_wind_turbine.f90'
!   include './g_move_wind_turbine-2.f90'
!   include './g_move_fish-sample.f90'
!   include './g_move_flapping.f90'

! ////////// Subroutine H //////////

   include './h_mpi_bcast_cal_data.f90'
   include './h_mpi_bcast_cal_update_data.f90'
   include './h_mpi_bcast_panel_data.f90'
   include './h_mpi_bcast_panel_update_data.f90'
   include './h_mpi_bcast_vor_data.f90'
   include './h_mpi_initial.f90'
   include './h_mpi_load_balance.f90'
   include './h_mpi_recv_vor_velo.f90'
   include './h_mpi_send_vor_velo.f90'

! ////////// Lapack lib. //////////

   include './LAPACK_library/dgemm.f'
   include './LAPACK_library/dger.f'
   include './LAPACK_library/dgesv.f'
   include './LAPACK_library/dgetf2.f'
   include './LAPACK_library/dgetrf.f'
   include './LAPACK_library/dgetrs.f'
   include './LAPACK_library/dlamch.f'
   include './LAPACK_library/dlaswp.f'
   include './LAPACK_library/dscal.f'
   include './LAPACK_library/dswap.f'
   include './LAPACK_library/dtrsm.f'
   include './LAPACK_library/idamax.f'
   include './LAPACK_library/ieeeck.f'
   include './LAPACK_library/ilaenv.f'
   include './LAPACK_library/iparmq.f'
   include './LAPACK_library/lsame.f'
   include './LAPACK_library/xerbla.f'


