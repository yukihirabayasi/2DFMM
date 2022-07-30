! ************************************************************************** 
!     program flow_field_contour
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING       
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  convert from s_ data to vtk format for paraview                           
!   (velo, pres & vorticity distribution at each time step can be displayed) 
!============================================================================

  include './includes/modules.f90'

  use global
  use parameter
  use math_constant
  implicit none
  include './includes/inc_2d'

   integer :: n_start
   integer :: n_end
   integer :: n_int

!   real(8) :: xi, xa, yi, ya
!   real(8) :: gid, gpx, gpy, gpz

!   character(50) :: nfile
!   character(40) :: fout_h
!   character(40) :: fout_t

   write(*,*) 'USE MIRROR IMAGE [0:no 1:yes]'
   read (*,*) i_mirror
   if ( i_mirror .eq. 1 ) then
     write(*,*) 'MIRRO DIRECTION'
     read (*,*) d_mirror
   end if
   write(*,*) 'START STEP NUMBER'
   read (*,*) n_start
   write(*,*) 'END STEP NUMBER'
   read (*,*) n_end
   write(*,*) 'SAVE INTERVAL'
   read (*,*) n_int

   open (11,file='./set_grid/grid_data.dat',status='old',form='formatted')
    read (11,*) nx
    read (11,*) ny
    read (11,*) xmin
    read (11,*) xmax
    read (11,*) ymin
    read (11,*) ymax
   close(11)

   do i = n_start, n_end, n_int

     call read_in (i)

     x_min = xmin + parts_position(1,1)
     x_max = xmax + parts_position(1,1)
     y_min = ymin + parts_position(2,1)
     y_max = ymax + parts_position(2,1)

     call get_velo_pres
     call data_out(i)
   end do

  stop
  end

  include './includes/d_biot.f90'
  include './includes/d_indus.f90'
  include './includes/d_biot_p.f90'


!************************************************************************** 
     subroutine read_in( int )
!************************************************************************** 

  use global
  use parameter
  use math_constant
  implicit none
  include './includes/inc_2d'

   integer :: int

   character(40) :: nfile
   character(40) :: fin_h
   character(40) :: fin_t

   fin_h = '../data_out/s_'
   fin_t = '_step.dat'

   write(nfile,'(a14,i5.5,a9)') fin_h, int, fin_t
   write(*,*) ''
   write(*,*) 'read_file   :   ', nfile

   open(10,file=nfile,status='old',form='unformatted')

    read(10) re, dt
    read(10) theta
    read(10) omz
    read(10) nt, time
    read(10) dh
    read(10) vn_vis
    read(10) ds_t
    read(10) dt_base
    read(10) accelerate
    read(10) velocity
    read(10) velocity_base
    read(10) uniform_velocity        ! add by K.F. (2012/08/13)
    read(10) uniform_velocity_base   ! add by K.F. (2012/08/13)
    read(10) u_accelerate            ! add by K.F. (2012/08/13)

    read(10) nw
    do i = 1, nw
      read(10) parts_position (1,i)
      read(10) parts_position (2,i)
      read(10) parts_position (3,i)
      read(10) pgw   (1,i)
      read(10) pgw   (2,i)
      read(10) pgw   (3,i)
      read(10) phase   (i)
      read(10) cir_r   (i)
      read(10) edge  (1,i)
      read(10) edge  (2,i)
    end do

    read(10) npanel
    do i = 1, npanel
      read(10) q        (i)
      read(10) ph       (i)
      read(10) ds       (i)
      read(10) cir_g    (i)
      read(10) uw     (1,i), uw     (2,i)
      read(10) poi  (1,1,i), poi  (2,1,i)
      read(10) poi  (1,2,i), poi  (2,2,i)
      read(10) poi_b(1,1,i), poi_b(2,1,i)
      read(10) poi_b(1,2,i), poi_b(2,2,i)
      read(10) poi_c  (1,i), poi_c  (2,i)
      read(10) poi_n  (1,i), poi_n  (2,i)
      read(10) vm   (1,1,i), vm   (2,1,i)
      read(10) vm   (1,2,i), vm   (2,2,i)
      read(10) am   (1,1,i), am   (2,1,i)
      read(10) am   (1,2,i), am   (2,2,i)
      read(10) panel_id (i)
      read(10) panel_n(1,i)
      read(10) panel_n(2,i)
    end do

    do i = 1, npanel
      read(10) nlayer   (i)
      read(10) lay_h    (i)
      read(10) ypls     (i)
      do j = 1, nlayer(i)
        read(10) lay_x  (1,j,i), lay_x  (2,j,i)
        read(10) lay_y  (1,j,i), lay_y  (2,j,i)
        read(10) lay_c  (1,j,i), lay_c  (2,j,i)
        read(10) lay_vx (1,j,i), lay_vx (2,j,i)
        read(10) lay_vy (1,j,i), lay_vy (2,j,i)
        read(10) lay_vc (1,j,i), lay_vc (2,j,i)
      end do
    end do

    read(10) nvor_b
    do i = 1, nvor_b
      read(10) vor_b    (1,i), vor_b    (2,i)
      read(10) vor_vb (1,1,i), vor_vb (2,1,i)
      read(10) vor_vb (1,2,i), vor_vb (2,2,i)
      read(10) vor_vb (1,3,i), vor_vb (2,3,i)
      read(10) cir_b      (i)
      read(10) cor_b      (i)
      read(10) blob_id    (i)
    end do

    read(10) nvor_s
    do i = 1, nvor_s
      read(10) vor_sc   (1,i), vor_sc   (2,i)
      read(10) vor_s  (1,1,i), vor_s  (2,1,i)
      read(10) vor_s  (1,2,i), vor_s  (2,2,i)
      read(10) vor_vs (1,1,i), vor_vs (2,1,i)
      read(10) vor_vs (1,2,i), vor_vs (2,2,i)
      read(10) vor_vs (1,3,i), vor_vs (2,3,i)
      read(10) cir_s      (i)
      read(10) cor_s      (i)
      read(10) vor_sr     (i)
      read(10) sheet_id   (i)
    end do

   close(10)

! //////////////  Sheet ---> Blob  //////////////

   do i = 1, nvor_s
     nvor_b = nvor_b + 1
     vor_b   (1,nvor_b) = vor_sc  (1,i)
     vor_b   (2,nvor_b) = vor_sc  (2,i)
     vor_vb(1,1,nvor_b) = vor_vs(1,1,i)
     vor_vb(1,2,nvor_b) = vor_vs(1,2,i)
     vor_vb(2,1,nvor_b) = vor_vs(2,1,i)
     vor_vb(2,2,nvor_b) = vor_vs(2,2,i)
     cir_b     (nvor_b) = cir_s     (i)
     cor_b     (nvor_b) = dsqrt(vor_sr(i)*2.0*cor_s(i)/pi)
     blob_id   (nvor_b) = sheet_id  (i)
   end do

  return
  end


! ***************************************************************************
     subroutine get_velo_pres
! ***************************************************************************
! -------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING       
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include './includes/inc_2d'

!  ////////// calculate velocity distribution ///////////

   uinf = uniform_velocity - velocity
   vinf = 0.0d0

   dx = dabs( x_max - x_min ) / dble( nx )
   dy = dabs( y_max - y_min ) / dble( ny )

   do i = 1, nx+3
   do j = 1, ny+3

     point(1,i,j) = x_min + dble( i-1 ) * dx
     point(2,i,j) = y_min + dble( j-1 ) * dy

     x1 = point(1,i,j)
     y1 = point(2,i,j)

     call cross1( x1,y1,icross )

     if ( icross .eq. 1 ) then
       velo(1,i,j) = 0.0d0
       velo(2,i,j) = 0.0d0
       pres  (i,j) = 0.0d0
       cycle
     end if

     call indus_s( x1,y1,u1,v1 )
     call indus_v( x1,y1,u2,v2 )
     call biot_b ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )

     call gh_b( x1,y1,1,bh_b )
     call gh_s( x1,y1,1,bh_s )
     call gh_r( x1,y1,1,bh_r )
     call gh_a( x1,y1,1,bh_a )

     vmx = uinf + ( -omz*y1 )
     vmy = vinf + (  omz*x1 )

     velo(1,i,j) =    vmx + ( u1 + u2 + u3 + u4 )
     velo(2,i,j) =    vmy + ( v1 + v2 + v3 + v4 )
     pres  (i,j) = -( bh_b + bh_s + bh_a + bh_r )

   end do
   end do

  return
  end


! ************************************************************************** 
      subroutine data_out( int )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include './includes/inc_2d'

   integer, parameter :: n4 = 4
   integer, parameter :: n9 = 9

   integer :: nc
   integer :: n_check
   integer :: int

   character(40) :: nfile
   character(40) :: fout_h
   character(40) :: fout_m
   character(40) :: fout_t

   dx = dabs( x_max - x_min ) / dble( nx )
   dy = dabs( y_max - y_min ) / dble( ny )

!  ////////// output data ///////////

   fout_h = '../data_out/body_'
   fout_m = '_'
   fout_t = '_step.vtk'

   do k = 1, nw

     write(nfile,'(A17,i2.2,A1,i5.5,A9)') fout_h, k, fout_m, int, fout_t

     open(11,file=nfile,status='unknown',form='formatted')

      write(11,'(A26)') '# vtk DataFile Version 2.0'
      write(11,'(A4)')  'Body'
      write(11,'(A5)')  'ASCII'
      write(11,'(A16)') 'DATASET POLYDATA'

      nc = 0

      do i = 1, npanel
        if ( panel_id(i) .eq. k ) nc = nc + 1
      end do

      write(11,'(A6,1x,i8,1x,A5)') 'POINTS', nc*4, 'float'

      do i = 1, npanel
        if ( panel_id(i) .eq. k ) then
          write(11,'(3(f16.8,2x))') poi(1,1,i), poi(2,1,i), 0.0d0
          write(11,'(3(f16.8,2x))') poi(1,2,i), poi(2,2,i), 0.0d0
          write(11,'(3(f16.8,2x))') poi(1,2,i), poi(2,2,i), 1.0d0
          write(11,'(3(f16.8,2x))') poi(1,1,i), poi(2,1,i), 1.0d0
        end if
      end do

      n_check = 0

      write(11,'(A8,1x,i8,1x,i8)') 'POLYGONS', nc, nc*5

      do i = 1, npanel
        if ( panel_id(i) .eq. k ) then
          write(11,'(5i8)') n4, n_check, n_check+1, n_check+2, n_check+3
          do j = 1, 4
            n_check = n_check + 1
          end do
        end if
      end do

     close(11)

   end do

   fout_h = '../data_out/flow_contour_'
   fout_t = '_step.vtk'

   write(nfile,'(A25,i5.5,A9)') fout_h, int, fout_t

   write(*,*) '  output_file : ', nfile
   write(*,*) '  nvor        : ', nvor_b

   open(12,file=nfile,status='unknown',form='formatted')

    write(12,'(A26)')            '# vtk DataFile Version 2.0'
    write(12,'(A21)')            'velocity_distribution'
    write(12,'(A5)')             'ASCII'
    write(12,'(A25)')            'DATASET UNSTRUCTURED_GRID'

    write(12,'(A1)')             ' '
    write(12,'(A6,2x,i8,2x,A5)') 'POINTS', nx*ny*4, 'float'

    do i = 1, nx
    do j = 1, ny
      write(12,'(3(e16.8,2x))')  point(1,i,  j  ), point(2,i,  j  ), 0.0d0
      write(12,'(3(e16.8,2x))')  point(1,i+1,j  ), point(2,i+1,j  ), 0.0d0
      write(12,'(3(e16.8,2x))')  point(1,i+1,j+1), point(2,i+1,j+1), 0.0d0
      write(12,'(3(e16.8,2x))')  point(1,i,  j+1), point(2,i,  j+1), 0.0d0
    end do
    end do

    n_check = 0

    write(12,'(A1)')             ' '
    write(12,'(A5,2x,i8,2x,i8)') 'CELLS', nx*ny, nx*ny*5

    do i = 1, nx*ny
      write(12,'(i1,2x,4(i8,2x))')  &
     &          n4, n_check, n_check+1, n_check+2, n_check+3
      do j = 1, 4
        n_check = n_check + 1
        end do
      end do

    write(12,'(A1)')             ' '
    write(12,'(A10,2x,i8)')      'CELL_TYPES', nx*ny

    do i = 1, nx*ny
      write(12,'(i1)') n9
    end do

    write(12,'(A1)')             ' '
    write(12,'(A10,2x,i8)')      'POINT_DATA', nx*ny*4
    write(12,'(A1)')             ' '
    write(12,'(A25)')            'VECTORS Velocity float'

    do i = 1, nx
    do j = 1, ny
      write(12,'(3(e16.8,2x))')  velo(1,i,  j  ), velo(2,i,  j  ), 0.0d0
      write(12,'(3(e16.8,2x))')  velo(1,i+1,j  ), velo(2,i+1,j  ), 0.0d0
      write(12,'(3(e16.8,2x))')  velo(1,i+1,j+1), velo(2,i+1,j+1), 0.0d0
      write(12,'(3(e16.8,2x))')  velo(1,i,  j+1), velo(2,i,  j+1), 0.0d0
    end do
    end do

    write(12,'(A25)')            'SCALARS Prresure float'
    write(12,'(A20)')            'LOOKUP_TABLE default'

    do i = 1, nx
    do j = 1, ny
      write(12,'(e16.8)') pres(i,  j  )
      write(12,'(e16.8)') pres(i+1,j  )
      write(12,'(e16.8)') pres(i+1,j+1)
      write(12,'(e16.8)') pres(i,  j+1)
    end do
    end do

    write(12,'(A25)')            'SCALARS Vorticity float'
    write(12,'(A20)')            'LOOKUP_TABLE default'

    do i = 2, nx+1
    do j = 2, ny+1
      write(12,'(e16.8)')  ( velo(2, i+1, j  ) - velo(2, i-1, j  ) ) / ( 2.0d0*dx ) &
     &                   - ( velo(1, i  , j+1) - velo(1, i  , j-1) ) / ( 2.0d0*dy )
      write(12,'(e16.8)')  ( velo(2, i+2, j  ) - velo(2, i  , j  ) ) / ( 2.0d0*dx ) &
     &                   - ( velo(1, i+1, j+1) - velo(1, i+1, j-1) ) / ( 2.0d0*dy )
      write(12,'(e16.8)')  ( velo(2, i+2, j+1) - velo(2, i  , j+1) ) / ( 2.0d0*dx ) &
     &                   - ( velo(1, i+1, j+2) - velo(1, i+1, j  ) ) / ( 2.0d0*dy )
      write(12,'(e16.8)')  ( velo(2, i+1, j+1) - velo(2, i-1, j+1) ) / ( 2.0d0*dx ) &
     &                   - ( velo(1, i  , j+2) - velo(1, i  , j  ) ) / ( 2.0d0*dy )
    end do
    end do

   close(12)

  return
  end


! ************************************************************************** 
     subroutine cross1( xpc,ypc,icross )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  use global
  use parameter
  use math_constant
  implicit none
  include './includes/inc_2d'

   icross = 0
   dmin   = 1.0d10

   do i = 1, npanel

     poicx = poi_c(1,i)
     poicy = poi_c(2,i)
     dmin2 = dsqrt( ( xpc-poicx )**2 + ( ypc-poicy )**2 )

     if ( dmin2 .lt. dmin ) then
       dmin = dmin2
       ic = i
     end if

   end do

   poinx = poi_n(1,ic)
   poiny = poi_n(2,ic)
   poicx = poi_c(1,ic)
   poicy = poi_c(2,ic)

   ax1 = poi(1,1,ic) - xpc
   ay1 = poi(2,1,ic) - ypc
   ax2 = poi(1,2,ic) - xpc
   ay2 = poi(2,2,ic) - ypc

   ss  = ax1*ay2 - ay1*ax2
   hei = dabs( ss/ds(ic) )
   height = hei

   rx1 = poi(1,1,ic) - poicx
   ry1 = poi(2,1,ic) - poicy
   r1  = dsqrt( rx1**2 + ry1**2 )
   fg1 = ( ax1*rx1 + ay1*ry1 ) / r1

   if ( fg1 .lt. 0.0d0 ) goto 100

   rx2 = poi(1,2,ic) - poicx
   ry2 = poi(2,2,ic) - poicy
   r2  = dsqrt( rx2**2 + ry2**2 )
   fg2 = ( ax2*rx2 + ay2*ry2 ) / r2

   if ( fg2 .lt. 0.0d0 ) goto 100

   vss = dsqrt( ( xpc - poicx )**2 + ( ypc - poicy )**2 )
   vin = poinx*( xpc - poicx ) + poiny*( ypc - poicy )

   if ( vin/vss .ge. 0.0d0 ) then
     if ( ( hei .gt. 0.0d0 ) .and. ( hei .le. dh ) ) icross = 1
     if (   hei .le. 0.0d0 )                         icross = 1
   else
     icross = 1
   end if

   100 continue

  return
  end


