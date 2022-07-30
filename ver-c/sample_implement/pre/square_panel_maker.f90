! ************************************************************************* 
!     program square_panel_maker
! ************************************************************************* 
! ------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING      
! AGREEMENT.                                                                
! ------------------------------------------------------------------------- 

 implicit none

  integer :: nw
  integer :: nx
  integer :: ny
  integer :: i

  double precision :: poi (2,1001)
  double precision :: dx
  double precision :: dy
  double precision :: rx
  double precision :: ry

  character(80) :: file1

   write(*,*) 'OUTPUT FILE NAME ?'
   read (*,'(a80)') file1
   write(*,*) 'HORIZONTAL    SIDE  LENGTH RX = ?'
   read (*,*) rx
   write(*,*) 'HORIZONTAL    PANEL NUMBER NX = ?'
   read (*,*) nx
   write(*,*) 'PERPENDICULAR SIDE  LENGTH RY = ?'
   read (*,*) ry
   write(*,*) 'PERPENDICULAR PANEL NUMBER NY = ?'
   read (*,*) ny

   dx = rx / nx
   dy = ry / ny

   do i = 1, nx+1
     poi(1,i)           =  0.5d0*rx - (i-1)*dx
     poi(2,i)           =  0.5d0*ry
   end do

   do i = 1, ny
     poi(1,nx+1+i)      = -0.5d0*rx
     poi(2,nx+1+i)      =  0.5d0*ry - i*dy
   end do

   do i = 1, nx
     poi(1,nx+ny+1+i)   = -0.5d0*rx + i*dx
     poi(2,nx+ny+1+i)   = -0.5d0*ry
   end do

   do i = 1, ny
     poi(1,2*nx+ny+1+i) =  0.5d0*rx
     poi(2,2*nx+ny+1+i) = -0.5d0*ry + i*dy
   end do

   nw = 1
   open (11,file=file1,status='unknown',form='formatted')

   write(11,*) nw
   write(11,*) 2*(nx+ny)+1

   do i = 1, 2*(nx+ny)+1
     write(11,'(2(f12.8,2x))') poi(1,i), poi(2,i)
   end do

   close(11)

   stop

   end

