! ************************************************************************* 
!     program elliptical_wing_maker
! ************************************************************************* 
! ------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING      
! AGREEMENT.                                                                
! ------------------------------------------------------------------------- 

 implicit none

  integer :: nw
  integer :: n1
  integer :: n2
  integer :: n
  integer :: nl
  integer :: i

  double precision :: poi (2,5001)
  double precision :: d1
  double precision :: d2
  double precision :: l_x
  double precision :: l_y
  double precision :: l
  double precision :: dl_x
  double precision :: dl_y
  double precision :: angle
  double precision :: d_angle
  double precision :: angle_s
  double precision :: angle_e
  double precision :: d_panel
  double precision :: d_panel_1
  double precision :: d_panel_2
  double precision :: pi

  character(80) :: file1

   write(*,*) 'OUTPUT FILE NAME ?'
   read (*,'(a80)') file1
   write(*,*) 'OUTER DIAMETER D1 = ?'
   read (*,*) d1
   write(*,*) 'INNER DIAMETER D2 = ?'
   read (*,*) d2
   write(*,*) 'PANEL LENGTH FOR SIDE EDGE = ?'
   read (*,*) d_panel_1
   write(*,*) 'MAXIMUM PANEL LENGTH = ?'
   read (*,*) d_panel_2
   write(*,*) 'ANGLE = ?'
   read (*,*) angle

   pi    = 4.0d0 * datan( 1.0d0 )
   angle = angle*pi/180.0d0

   n1 = int( d1*angle / 2.0d0 / ( d_panel_1 + d_panel_2 ) )
   n2 = int( d2*angle / 2.0d0 / ( d_panel_1 + d_panel_2 ) )

   n  = 0

   angle_s = 0.5d0*( pi - angle )
   angle_e = 0.5d0*( pi + angle )

   n = n + 1

   poi(1,n) = 0.5d0*d1*dcos( angle_s )
   poi(2,n) = 0.5d0*d1*dsin( angle_s )

   d_angle = angle_s
   d_panel = 2.0d0 * ( d1*angle / 4.0d0 / dble( n1 ) - d_panel_1 )

   do i = 1, n1
     n = n + 1
     d_angle  = d_angle + 2.0d0 / d1 * ( d_panel_1 + d_panel*dble( i-1 )/dble( n1-1 ) )
     poi(1,n) = 0.5d0*d1*dcos( d_angle )
     poi(2,n) = 0.5d0*d1*dsin( d_angle )
     if ( i .eq. n1 ) write(*,*) d_angle - 0.5d0*pi
   end do

   do i = 1, n1
     n = n + 1
     d_angle  = d_angle + 2.0d0 / d1 * ( d_panel_1 + d_panel*dble( n1-i )/dble( n1-1 ) )
     poi(1,n) = 0.5d0*d1*dcos( d_angle )
     poi(2,n) = 0.5d0*d1*dsin( d_angle )
     if ( i .eq. n1 ) write(*,*) d_angle - angle_e
   end do

   l_x  = 0.5d0*d2*dcos( angle_e ) - 0.5d0*d1*dcos( angle_e )
   l_y  = 0.5d0*d2*dsin( angle_e ) - 0.5d0*d1*dsin( angle_e )
   l    = dsqrt( l_x**2 + l_y**2 )

   nl = int( l/d_panel_1 )

   dl_x = l_x / dble( nl )
   dl_y = l_y / dble( nl )

   do i = 1, nl
     n = n + 1
     poi(1,n) = 0.5d0*d1*dcos( angle_e ) + dl_x*dble( i )
     poi(2,n) = 0.5d0*d1*dsin( angle_e ) + dl_y*dble( i )
   end do

   d_angle = angle_e
   d_panel = 2.0d0 * ( d2*angle / 4.0d0 / dble( n2 ) - d_panel_1 )

   do i = 1, n2
     n = n + 1
     d_angle  = d_angle - 2.0d0 / d2 * ( d_panel_1 + d_panel*dble( i-1 )/dble( n2-1 ) )
     poi(1,n) = 0.5d0*d2*dcos( d_angle )
     poi(2,n) = 0.5d0*d2*dsin( d_angle )
     if ( i .eq. n2 ) write(*,*) d_angle - 0.5d0*pi
   end do

   do i = 1, n2
     n = n + 1
     d_angle  = d_angle - 2.0d0 / d2 * ( d_panel_1 + d_panel*dble( n2-i )/dble( n2-1 ) )
     poi(1,n) = 0.5d0*d2*dcos( d_angle )
     poi(2,n) = 0.5d0*d2*dsin( d_angle )
     if ( i .eq. n2 ) write(*,*) d_angle - angle_s
   end do

   l_x  = 0.5d0*d1*dcos( angle_s ) - 0.5d0*d2*dcos( angle_s )
   l_y  = 0.5d0*d1*dsin( angle_s ) - 0.5d0*d2*dsin( angle_s )
   l    = dsqrt( l_x**2 + l_y**2 )

   nl = int( l/d_panel_1 )

   dl_x = l_x / dble( nl )
   dl_y = l_y / dble( nl )

   do i = 1, nl
     n = n + 1
     poi(1,n) = 0.5d0*d2*dcos( angle_s ) + dl_x*dble( i )
     poi(2,n) = 0.5d0*d2*dsin( angle_s ) + dl_y*dble( i )
   end do

   do i = 1, n
     poi(1,i) = poi(1,i) - 0.5d0*d1*dcos( angle_s )
     poi(2,i) = poi(2,i) - 0.5d0*d1*dsin( angle_s )
   end do

   do i = 1, n
     poi(1,i) = poi(1,i) / ( d1*pi*angle/2.0d0/pi )
     poi(2,i) = poi(2,i) / ( d1*pi*angle/2.0d0/pi )
   end do

   nw = 1
   open (11,file=file1,status='unknown',form='formatted')

   write(11,*) nw
   write(11,*) n-1

   do i = 1, n
     write(11,'(2(f12.8,2x))') poi(1,i), poi(2,i)
   end do

   close(11)

   stop

   end

