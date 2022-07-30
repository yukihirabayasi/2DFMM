! ************************************************************************* 
!     program naca_wing_panel_maker
! ************************************************************************* 
! ------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING      
! AGREEMENT.                                                                
! ------------------------------------------------------------------------- 

 implicit none

  integer :: nw
  integer :: n
  integer :: nc
  integer :: m
  integer :: i

  double precision :: upper(2,1001)
  double precision :: lower(2,1001)
  double precision :: t
  double precision :: x
  double precision :: x1
  double precision :: y1
  double precision :: mx
  double precision :: mc
  double precision :: yt
  double precision :: ym
  double precision :: dy_dx
  double precision :: x0
  double precision :: y0
  double precision :: dl
  double precision :: pi
  double precision :: theta
  double precision :: dtheta
  double precision :: angle

  character(80) :: file1

    go to 100

110 write(*,*) ' '
    write(*,*) 'To EXIT is [Ctrl] + [C]'
    write(*,*) ' '

100 write(*,*) ' '

    write(*,*) '   ex) NACA2412   '
    write(*,*) '   M=[2], P=[4], T=[12]'
    write(*,*) ' '
    write(*,*) 'OUTPUT FILE NAME ?'
    read (*,'(a80)') file1
    write(*,*) 'PANEL NUMBER = ?'
    read (*,*) n
    write(*,*) 'MAX CAMBER NUMBER [M]= ?'
    read (*,*) mc
    write(*,*) 'X-COODINATE OF MAX CAMBER NUMBER [P] = ?'
    read (*,*) mx
    write(*,*) 'WING THICKNESS [T] = ?'
    read (*,*) t

   pi = 4.0d0 * datan( 1.0d0 )

   mc     = mc * 0.01d0
   t      = t  * 0.01d0
   mx     = mx * 0.1d0
   dtheta = 2.0d0*pi/n

   m  = n / 2
   nc = m + 1

   do i = 1, m+1
     if ( i .le. m/2 ) then
       x = 1.0d0 - ( i-1 ) * 1.0d0 / m
     else
       theta = ( i-1 ) * dtheta
       x = 0.5d0 + 0.5d0 * dcos( theta )
     end if

     yt = t / 0.2d0 *                    &
   &              (  0.29690d0*dsqrt(x)  &
   &               - 0.12600d0*x**1      &
   &               - 0.35160d0*x**2      &
   &               + 0.28430d0*x**3      &
   &               - 0.10150d0*x**4 )

     if( x .lt. mx ) then
       ym    = mc * ( x/( mx**2 ) ) * ( 2.0d0*mx - x )
       dy_dx = 2.0d0 * mc *( (mx -x)/(mx**2) )
     else
       ym    = mc * ( (1.0d0 - x)/( (1.0d0 - mx)**2 ) ) * ( 1.0d0 + x - 2.0d0*mx )
       dy_dx = ( 2.0d0 * mc *(mx -x) )/( 1.0d0 - mx**2 )
     end if

     upper(1,i)    = x  - yt * dsin( datan( dy_dx ) )
     upper(2,i)    = ym + yt * dcos( datan( dy_dx ) )
     lower(1,nc-i) = x  + yt * dsin( datan( dy_dx ) )
     lower(2,nc-i) = ym - yt * dcos( datan( dy_dx ) )

     if ( i .eq. 1 ) then
       upper(2,i)    = 0.0d0
       lower(2,nc-i) = 0.0d0
     end if

   end do

   write(*,*) '[x0, y0] = '
   read (*,*)   x0, y0

   do i = 1, m+1
     upper(1,i) = upper(1,i) + x0
     upper(2,i) = upper(2,i) + y0
     lower(1,i) = lower(1,i) + x0
     lower(2,i) = lower(2,i) + y0
   end do

   write(*,*) 'ATTACK ANGLE [DEG]'
   read (*,*) angle

   angle = angle*pi/180.0d0

   do i = 1, m+1
     x1 = upper(1,i)
     y1 = upper(2,i)
     upper(1,i) =  x1*dcos( angle ) + y1*dsin( angle )
     upper(2,i) = -x1*dsin( angle ) + y1*dcos( angle )
   end do

   do i = 1, m+1
     x1 = lower(1,i)
     y1 = lower(2,i)
     lower(1,i) =  x1*dcos( angle ) + y1*dsin( angle )
     lower(2,i) = -x1*dsin( angle ) + y1*dcos( angle )
   end do

   write(*,*) 'CHORD LENGTH [L]'
   read (*,*) dl

   do i = 1, m+1
     upper(1,i) = dl * upper(1,i)
     upper(2,i) = dl * upper(2,i)
     lower(1,i) = dl * lower(1,i)
     lower(2,i) = dl * lower(2,i)
   end do

   nw = 1
   open (11,file=file1,status='unknown',form='formatted')

   write(11,*) nw
   write(11,*) m*2

   do i = 1, m+1
     write(11,'(2(f12.8,2x))') upper(1,i), upper(2,i)
   end do

   do i = 1, m
     write(11,'(2(f12.8,2x))') lower(1,i), lower(2,i)
   end do

   close(11)

   go to 110

   stop

   end

