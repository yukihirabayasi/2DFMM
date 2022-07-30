! ************************************************************************* 
!     program circle_panel_maker
! ************************************************************************* 
! ------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING      
! AGREEMENT.                                                                
! ------------------------------------------------------------------------- 

 implicit none

  integer :: i
  integer :: nw, n_panel

  double precision :: poi(2,1001)
  double precision :: theta, pi
  double precision :: ra, rb

  character(80) :: nfile

   pi = 4.0d0 * datan( 1.0d0 )
   nw = 1

   write(*,*) 'OUTPUT FILE NAME ?'
   read (*,*) nfile
   write(*,*) 'PANEL NUMBER = ?'
   read (*,*) n_panel
   write(*,*) 'RADIUS 1 (HORIZONTAL)    = ?'
   read (*,*) ra
   write(*,*) 'RADIUS 2 (PERPENDICULAR) = ?'
   read (*,*) rb

   open(10,file=nfile,status='unknown',form='formatted')

   write(10,*) nw
   write(10,*) n_panel

   do i = 1, n_panel+1
     theta = 2.0d0 * pi * ( dble( i )-1.0d0 ) / dble( n_panel )
     poi(1,i) = ra * dcos( theta )
     poi(2,i) = rb * dsin( theta )
     write(10,'(2(f12.8,2x))') poi(1,i), poi(2,i)
   end do

   close(10)

   stop

   end

