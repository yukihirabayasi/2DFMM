! ************************************************************************** 
!  convert from s_ data to vtk format for paraview                           
!   (vortex elements distribution at each time step can be displayed)        
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING       
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  include './includes/modules.f90'

  use global
  use parameter
  use math_constant
  implicit none

   integer :: n_start
   integer :: n_end
   integer :: n_int
   integer :: i

   write(*,*) 'START STEP NUMBER'
   read (*,*) n_start
   write(*,*) 'END STEP NUMBER'
   read (*,*) n_end
   write(*,*) 'SAVE INTERVAL'
   read (*,*) n_int

   do i = n_start, n_end, n_int
     call read_in  (i)
     call data_out (i)
   end do

  stop
  end


! ************************************************************************** 
     subroutine read_in (int)
! ************************************************************************** 

  use global
  use parameter
  use math_constant
  implicit none

   character discrption*80

   character(40) :: nfile
   character(40) :: fin_h
   character(40) :: fin_t
   character(40) :: body_file

   integer :: int
   integer :: i, j

   fin_h = '../data_out/s_'
   fin_t = '_step.dat'

   write(nfile,'(a14,i5.5,a9)') fin_h,int,fin_t
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


! ************************************************************************** 
     subroutine data_out (int)
! ************************************************************************** 

  use global
  use parameter
  use math_constant
  implicit none

   character(40) :: nfile
   character(40) :: fout_h
   character(40) :: fout_m
   character(40) :: fout_t
   integer :: n2
   integer :: n4
   integer :: nc
   integer :: n_check
   integer :: i, j, k
   integer :: int

   fout_h = '../data_out/body_'
   fout_m = '_'
   fout_t = '_step.vtk'

   do k = 1, nw

     write(nfile,'(a17,i2.2,a1,i5.5,a9)')fout_h,k,fout_m,int,fout_t

     open  (11,file=nfile,status='unknown',form='formatted')
      write(11,'(a26)') '# vtk DataFile Version 2.0' 
      write(11,'(a4)')  'Body'
      write(11,'(a5)')  'ASCII'
      write(11,'(a16)') 'DATASET POLYDATA'

      nc = 0
      do i = 1, npanel
        if ( panel_id(i) .eq. k ) nc = nc + 1
      end do

      write(11,'(a6,x,i8,x,a5)') 'POINTS',nc*4,'float'

      do i = 1, npanel
        if ( panel_id(i) .eq. k ) then
          write(11,'(f16.8,2x,f16.8,2x,f16.8)') poi(1,1,i), poi(2,1,i), 0.0d0
          write(11,'(f16.8,2x,f16.8,2x,f16.8)') poi(1,2,i), poi(2,2,i), 0.0d0
          write(11,'(f16.8,2x,f16.8,2x,f16.8)') poi(1,2,i), poi(2,2,i), 1.0d0
          write(11,'(f16.8,2x,f16.8,2x,f16.8)') poi(1,1,i), poi(2,1,i), 1.0d0
        end if
      end do

      n4      = 4
      n_check = 0
      write(11,'(a8,x,i8,x,i8)') 'POLYGONS', nc, nc*5

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

   fout_h = '../data_out/vorton_'
   fout_t = '_step.vtk'
   write(nfile,'(a19,i5.5,a9)')fout_h,int,fout_t

   write(*,*) 'output_file : ',nfile
   write(*,*) 'nvor        : ',nvor_b

   open (12,file=nfile,status='unknown',form='formatted')

    write(12,'(a26)') '# vtk DataFile Version 2.0'
    write(12,'(a27)') 'vortex_element_distribution'
    write(12,'(a5)')  'ASCII'

    write(12,'(a1)')  ''
    write(12,'(a16)') 'DATASET POLYDATA'
    write(12,'(a1)')  ''
    write(12,'(a6, 2x, i8, 2x, a5)') 'POINTS',nvor_b*2,'float'

    do i = 1, nvor_b
      write(12,'(e16.8, 2x, e16.8, 2x, e16.8)') vor_b(1,i), vor_b(2,i), 0.0d0
      write(12,'(e16.8, 2x, e16.8, 2x, e16.8)') vor_b(1,i), vor_b(2,i), 1.0d0
    end do

    write(12,'(a1)') ''
    write(12,'(a5, 2x, i8, 2x, i8)') 'LINES', nvor_b, nvor_b*3

    n2 = 2
    do i = 1, nvor_b
      write(12,'(i8, 2x, i8, 2x, i8)') n2, (i-1)*2, (i-1)*2+1
    end do

    write(12,'(a1)')       ''
    write(12,'(a10, 2x, i8)') 'POINT_DATA',nvor_b*2

    write(12,'(a1)')       ''
    write(12,'(a21)')      'SCALARS omega float'
    write(12,'(a20)')      'LOOKUP_TABLE default' 

    do i = 1, nvor_b
      write(12,'(e16.8)') cir_b(i)/pi/cor_b(i)**2
      write(12,'(e16.8)') cir_b(i)/pi/cor_b(i)**2
    end do

    write(12,'(a1)')       ''
    write(12,'(a21)')      'SCALARS ID float'
    write(12,'(a20)')      'LOOKUP_TABLE default' 

    do i = 1, nvor_b
      write(12,'(i3)') blob_id(i)
      write(12,'(i3)') blob_id(i)
    end do

    write(12,'(A1)')                 ' '
    write(12,'(A25)')                'VECTORS velocity float'

    do i = 1, nvor_b
      write(12,'(e16.8)')            vor_vb(1,1,i), vor_vb(2,1,i), 0.0d0
      write(12,'(e16.8)')            vor_vb(1,1,i), vor_vb(2,1,i), 0.0d0
    end do

   close(12)

  return
  end

