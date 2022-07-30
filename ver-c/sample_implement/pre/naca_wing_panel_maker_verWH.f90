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

  double precision,allocatable :: upper(:,:)
  double precision,allocatable :: lower(:,:)
  double precision :: x0, y0
  double precision :: x1, y1
  double precision :: mc, mp, mt
  double precision :: xc
  double precision :: yc, yt
  double precision :: dy_dx
  double precision :: pi
  double precision :: theta
  double precision :: dtheta
  double precision :: angle

  character(80) :: nfile

   print *, ' '
   print *, '  ex) NACA2412'
   print *, '  M=[2], P=[4], T=[12]'
   print *, ' '
   print *, 'OUTPUT FILE NAME ?'
   read  '(A)', nfile
   print *, 'PANEL NUMBER = ?'
   read  *, n
   print *, 'MAX CAMBER NUMBER [M]= ?'
   read  *, mc
   print *, 'X-COODINATE OF MAX CAMBER NUMBER [P] = ?'
   read  *, mp
   print *, 'WING THICKNESS [T] = ?'
   read  *, mt

   mc = mc * 0.01d0
   mp = mp * 0.10d0
   mt = mt * 0.01d0

   pi = 4.0d0 * datan( 1.0d0 )
   dtheta = pi / n

   m  = n / 2
   nc = m + 1

   allocate( upper(2,nc), lower(2,nc) )

   do i = 1, nc

!  ////////// Set X-Coodinate ///////////

     theta = 0.5d0*pi + dble( i-1 )*dtheta
     xc    = 1.0d0 + dcos( theta )

!  ////////// Cal. Thickness distribution ///////////

     yt = mt / 0.2d0 * (  0.29690d0 * dsqrt( xc )  &
   &                    - 0.12600d0 * xc           &
   &                    - 0.35160d0 * xc**2.0d0    &
   &                    + 0.28430d0 * xc**3.0d0    &
   &                    - 0.10360d0 * xc**4.0d0 )	! closed traling edge
!   &                    - 0.10150d0 * xc**4.0d0 )	! opened traling edge

!  ////////// Cal. Camber and Gradient  ///////////

     if ( xc .lt. mp ) then
       yc    = mc/(mp**2) * ( 2.0d0*mp*xc - xc**2 )
       dy_dx = 2.0d0*mc/(mp**2) * ( mp - xc )
     else
       yc    = mc/( (1.0d0-mp)**2 ) * ( 1.0d0 - 2.0d0*mp + 2.0d0*mp*xc - xc**2 )
       dy_dx = 2.0d0*mc /( (1.0d0-mp)**2 ) * ( mp - xc )
     end if

!  ////////// Cal. Upper & Lower surface point  ///////////

     upper(1,i)      = xc - yt * dsin( datan( dy_dx ) )
     upper(2,i)      = yc + yt * dcos( datan( dy_dx ) )
     lower(1,nc+1-i) = xc + yt * dsin( datan( dy_dx ) )
     lower(2,nc+1-i) = yc - yt * dcos( datan( dy_dx ) )

     if ( i .eq. 1 ) then
       upper(2,1)  = 0.0d0
       lower(2,nc) = 0.0d0
     end if

   end do

   print *, '[x0, y0] = '
   read *,   x0,  y0

   do i = 1, nc
     upper(1,i) = upper(1,i) + x0
     upper(2,i) = upper(2,i) + y0
     lower(1,i) = lower(1,i) + x0
     lower(2,i) = lower(2,i) + y0
   end do

   print *, 'ATTACK ANGLE [DEG]'
   read *, angle

   angle = angle*pi/180.0d0

   do i = 1, nc
     x1 = upper(1,i)
     y1 = upper(2,i)
     upper(1,i) =  x1*dcos( angle ) + y1*dsin( angle )
     upper(2,i) = -x1*dsin( angle ) + y1*dcos( angle )
     x1 = lower(1,i)
     y1 = lower(2,i)
     lower(1,i) =  x1*dcos( angle ) + y1*dsin( angle )
     lower(2,i) = -x1*dsin( angle ) + y1*dcos( angle )
   end do

   nw = 1

   open (11,file=nfile,status='unknown',form='formatted')

    write(11,*) nw
    write(11,*) n

    do i = 1, m
      write(11,'(2(f12.8,2x))') upper(1,i), upper(2,i)
    end do

    do i = 1, nc
      write(11,'(2(f12.8,2x))') lower(1,i), lower(2,i)
    end do

   close(11)

   deallocate( upper, lower )

  stop
  end


