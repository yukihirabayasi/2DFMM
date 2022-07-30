! ************************************************************************** 
     subroutine initial_check
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

   character(8)  :: date
   character(4)  :: date1
   character(2)  :: date2
   character(2)  :: date3
   character(10) :: time1
   character(26) :: ndate1
   character(26) :: ndate2
   character(26) :: fdate1
   character(26) :: fdate2

   write(*,*) ' '
   write(*,*) '*** Code License check ***'

   call date_and_time( date,time1 )

   ndate1 = date(1:4)//date(5:6)//date(7:8)
   date1  = '2031'
   date2  = '03'
   date3  = '31'
   ndate2 = date1//date2//date3

   fdate1 = ' System Date :  '//date(1:4)//'/'//date(5:6)//'/'//date(7:8)
   fdate2 = ' Valid Until :  '//date1    //'/'//date2    //'/'//date3

   write(*,*) fdate1
   write(*,*) fdate2

   if ( ndate1 .ge. ndate2 ) then
     write(*,*) '<<< This code is already expired !! >>>'
     stop
   end if

   open (12,file='../data_out/check_job_state.dat',status='unknown')
    write(12,*)                    ''
    write(12,*)                    '*=========================================*'
    write(12,*)                    '*  Numerical Simulation of Unsteady Flow  *'
    write(12,*)                    '*  around bodies by a 2-D Vortex Method   *'
    write(12,*)                    '*  for Parallel Processor with MPI        *'
    write(12,*)                    '*    coded by Kota Fukuda (2004/03/14)    *'
    write(12,*)                    '*=========================================*'
    write(12,*)                    ''
    write(12,*)                    ''
    write(12,*)                    '**** Code License check ****'
    write(12,*)                    fdate1
    write(12,*)                    fdate2
    write(12,*)                    ''
    write(12,*)                    ''
    write(12,*)                    '**** start time ****'
    write(12,*)                    fdate1
    write(12,*)                    ' System time :   ' &
   &                               //time1(1:2)//':'  &
   &                               //time1(3:4)//'-'  &
   &                               //time1(5:6),'s'
    write(12,*)                    ''
    write(12,*)                    ''
    write(12,*)                    '*** Calculation condition ***'
    write(12, '(1x, A, f10.2)')    ' Reynolds number   = ', re
    write(12, '(1x, A, f10.5)')    ' Time interval     = ', dt
    write(12, '(1x, A, f10.5)')    ' Attack angle      = ', theta
    write(12, '(1x, A, f10.5)')    ' Calculation Space = ', cut_r
    write(12, '(1x, A, f10.5)')    ' Uniform flow (x)  = ', uinf
    write(12, '(1x, A, f10.5)')    ' Uniform flow (y)  = ', vinf
    write(12, '(1x, A, f10.5)')    ' Angular velocity  = ', omz
    write(12, '(1x, A, f10.5)')    ' dh                = ',dh
    write(12, '(1x, A, f10.5)')    ' rotation radius   =', r_rot
    write(12, '(1x, A, f10.5)')    ' rotation velosity =', v_rot
    write(12, '(1x, A, f10.5)')    ' tip speed ratio   =', tsr
    write(12, '(1x, A, f10.5)')    ' cycle_time        =', cycle_time
    write(12, '(1x, A, f10.5)')    ' pitch_angle_up    =', pitch_angle_base_up
    write(12, '(1x, A, f10.5)')    ' pitch_angle_down  =', pitch_angle_base_down
    write(12, '(1x, A, f10.5)')    ' interval for heaving           =',dt_t
    write(12, '(1x, A, f10.5)')    ' interval for pitching   (dt_r) =',dt_r
    write(12, '(1x, A, f10.5)')    ' interval for pitching (dt_lag) =',dt_lag

    if ( i_continue .eq. 1 ) write(12,*) ' continued file    = ', file2
    write(12,*)                    ''
    write(12,*)                    ''
    write(12,*)                    '*** Panel data ***'
    write(12,*)                    ' Panel data   : ',file1
    write(12, '(2x, A, i2)')       'Parts number : ', nw
    do i = 1, nw 
      n = 0
      do j = 1, npanel
        if( panel_id(j) .eq. i ) n = n + 1
      end do
      write(12,'(2x,A,i2,A,i5,A)')   '   ID : ',i,'   includes ',n,' elements'
    end do
    write(12,*)                    ''
    write(12,*)                    ''
   close(12)

  return
  end


