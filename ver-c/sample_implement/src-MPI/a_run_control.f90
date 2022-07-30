! ************************************************************************** 
     subroutine create_run_check_file
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  integer  :: n_check_run

   n_check_run = 0

   open (18,file='../data_out/check_job_control.dat',status='unknown',form='formatted')
    write(18,'(i5,10x,A)') n_check_run,'/RUN CONTROLER [O:RUN / 1:STOP] '
   close(18)

  return
  end


! ************************************************************************** 
     subroutine run_check_and_close_run
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

   character (8)  :: date1
   character(10)  :: time1
   integer  :: n_check_run

   open (18,file='../data_out/check_job_control.dat',status='unknown',form='formatted')
    read (18,*) n_check_run
   close(18)

!=================================
   if ( n_check_run .eq. 1 ) then 
!=================================

   write(*, '(A)') ''
   write(*, '(A)') ''
   write(*, '(A)') ''
   write(*, '(A)') '  <<< This process is being finished.  Please wait !! >>>'

   call file_close

   write(*, '(A)') ''
   write(*, '(A)') '  ***** The job was successfully stopped !! *****'
   write(*, '(A)') ''
   write(*, '(A)') ''

   call date_and_time( date1,time1 )

   open (19,file='../data_out/check_job_state.dat' &
  &               ,position='append',status='unknown',form='formatted')
    write(19,*)           ''
    write(19,'(A)')       ' ***** The job was successfully stoped *****'
    write(19,*)           ' System Date :  ' &
   &                       //date1(1:4)//'/'//date1(5:6)//'/'//date1(7:8)
    write(19,*)           ' System time :   ' &
   &                       //time1(1:2)//':'//time1(3:4)//'-'//time1(5:6),'s'
    write(19,'(A,i10)')    '  time step   :  ',nt - 1
    write(19,'(A,f10.6)')  '  time        :  ',time - dt
   close(19)
   stop

!=========
   end if 
!=========

  return
  end


