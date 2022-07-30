! ************************************************************************* 
!     program panel_move
! ************************************************************************* 
! ------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING      
! AGREEMENT.                                                                
! ------------------------------------------------------------------------- 

 implicit none

  integer :: i
  integer :: nw, n

  double precision :: poi(2,1001)
  double precision :: pi
  double precision :: x0, y0
  double precision :: x1, y1
  double precision :: theta
  double precision :: dtheta

  character(80) :: file1
  character(80) :: file2

   pi = 4.0d0 * datan( 1.0d0 )

   write(*,*) 'INPUT FILE NAME ?'
   read (*,'(a80)') file1

   open (14,file=file1,status='old',form='formatted')

   read(14,*) nw
   read(14,*) n

   do i = 1, n+1
     read(14,*) poi(1,i), poi(2,i)
   end do

   close(14)

   write(*,*) 'O-POINT[X0,Y0] = ?'
   read (*,*) x0, y0
   write(*,*) 'ANGLE [DEG] = ?'
   read (*,*) theta

   dtheta = 2.0d0*pi*theta/360.0d0

   do i = 1, n+1
     x1 = poi(1,i)
     y1 = poi(2,i)
     poi(1,i) = x1*dcos( dtheta ) - y1*dsin( dtheta ) + x0
     poi(2,i) = x1*dsin( dtheta ) + y1*dcos( dtheta ) + y0
   end do

   write(*,*) 'OUTPUT FILE NAME ?'
   read(*,'(a80)') file2

   open (15,file=file2,status='unknown',form='formatted')

   write(15,*) nw
   write(15,*) n

   do i = 1, n+1
     write(15,'(2(f12.8,2x))') poi(1,i), poi(2,i)
   end do

   close(15)

   stop

   end

