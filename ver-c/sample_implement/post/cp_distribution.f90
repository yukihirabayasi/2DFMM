! ************************************************************************** 
!    program Cp distribution
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING       
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Read p_file and calculate non-dimensional presure                         
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

   character(40) nfile
   character(40) fout_h
   character(40) fout_t

   write(*,*) 'START STEP NUMBER'
   read (*,*) n_start
   write(*,*) 'END STEP NUMBER'
   read (*,*) n_end
   write(*,*) 'SAVE INTERVAL'
   read (*,*) n_int

   do i = n_start, n_end, n_int
     call pos_read(i)
     call read_in (i)
     call data_out(i)
   end do

  stop
  end


! ************************************************************************** 
     subroutine read_in( int )
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

   integer :: int

   character(40) :: nfile
   character(40) :: fin_h
   character(40) :: fin_t
   character(40) :: fout_h
   character(40) :: fout_t

   fin_h = '../data_out/s_'
   fin_t = '_step.dat'

   write(nfile,'(a14,i5.5,a9)') fin_h, int, fin_t
   write(*,*)  '  read_file   : ', nfile

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

   fout_h = '../data_out/p_'
   fout_t = '_step.dat'
   write(nfile,'(a14,i5.5,a9)') fout_h, int, fout_t

   open(10,file=nfile,status='old',form='unformatted')
    read(10) npanel
    do i = 1, npanel
      read(10) poi_c(1,i), poi_c(2,i), cp_t(i), cf(i)
    end do
   close(10)

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

   integer :: int

   character(40) :: nfile
   character(40) :: fout_h
   character(40) :: fout_t

   fout_h = '../data_out/cp_dis_'
   fout_t = '_step.dat'

   write(nfile,'(a19,i5.5,a9)') fout_h, int, fout_t
   write(*,*) '  output file : ', nfile

   open(20,file=nfile,status='unknown',form='formatted')
    write(20,*) '   poi_c_x,         non-cp,         cp'
    do i = 1, npanel
      write(20,'(3f18.12)') poi_c(1,i)-gpx, cp_t(i)/velocity**2, cp_t(i)
    end do
   close(20)

  return
  end


! ************************************************************************** 
     subroutine pos_read( int )
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

   integer :: int

   character(40) :: fout_h
   character(40) :: fout_t
   character(50) :: nfile

   fout_h = '../data_out/panel_position_'
   fout_t = '_step.dat'

   write(nfile,'(a27,i5.5,a9)') fout_h,int,fout_t
   write(*,*)  ' '
   write(*,*)  '  read file   : ', nfile

   open (11,file=nfile,status='old',form='formatted')
    read (11,*) gid, gpx, gpy, gpz
   close(11)

  return
  end


