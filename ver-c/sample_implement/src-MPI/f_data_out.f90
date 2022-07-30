! ************************************************************************** 
     subroutine data_out
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

   character(40) :: nfile
   character(40) :: fout_h
   character(40) :: fout_t

   fout_h = '../data_out/s_'
   fout_t = '_step.dat'

   write(nfile,'(a14,i5.5,a9)') fout_h, nt, fout_t

   write(*,*) ''
   write(*,*) '  data output : ', nfile

   nout = 16
   open(nout,file=nfile,status='unknown',form='unformatted')

    write(nout) re, dt
    write(nout) theta
    write(nout) omz
    write(nout) nt, time
    write(nout) dh
    write(nout) vn_vis
    write(nout) ds_t
    write(nout) dt_base
    write(nout) accelerate
    write(nout) velocity
    write(nout) velocity_base
    write(nout) uniform_velocity        ! add by K.F. (2012/08/13)
    write(nout) uniform_velocity_base   ! add by K.F. (2012/08/13)
    write(nout) u_accelerate            ! add by K.F. (2012/08/13)

    write(nout) nw

    do i = 1, nw
      write(nout) parts_position (1,i)
      write(nout) parts_position (2,i)
      write(nout) parts_position (3,i)
      write(nout) pgw   (1,i)
      write(nout) pgw   (2,i)
      write(nout) pgw   (3,i)
      write(nout) phase   (i)
      write(nout) cir_r   (i)
      write(nout) edge  (1,i)
      write(nout) edge  (2,i)
    end do
 
    write(nout) npanel

    do i = 1 , npanel
      write(nout) q        (i)
      write(nout) ph       (i)
      write(nout) ds       (i)
      write(nout) cir_g    (i)
      write(nout) uw     (1,i), uw     (2,i)
      write(nout) poi  (1,1,i), poi  (2,1,i)
      write(nout) poi  (1,2,i), poi  (2,2,i)
      write(nout) poi_b(1,1,i), poi_b(2,1,i)
      write(nout) poi_b(1,2,i), poi_b(2,2,i)
      write(nout) poi_c  (1,i), poi_c  (2,i)
      write(nout) poi_n  (1,i), poi_n  (2,i)
      write(nout) vm   (1,1,i), vm   (2,1,i)
      write(nout) vm   (1,2,i), vm   (2,2,i)
      write(nout) am   (1,1,i), am   (2,1,i)
      write(nout) am   (1,2,i), am   (2,2,i)
      write(nout) panel_id (i)
      write(nout) panel_n(1,i)
      write(nout) panel_n(2,i)
    end do

    do i = 1, npanel
      write(nout) nlayer   (i)
      write(nout) lay_h    (i)
      write(nout) ypls     (i)
      do j = 1, nlayer(i)
        write(nout) lay_x  (1,j,i), lay_x  (2,j,i)
        write(nout) lay_y  (1,j,i), lay_y  (2,j,i)
        write(nout) lay_c  (1,j,i), lay_c  (2,j,i)
        write(nout) lay_vx (1,j,i), lay_vx (2,j,i)
        write(nout) lay_vy (1,j,i), lay_vy (2,j,i)
        write(nout) lay_vc (1,j,i), lay_vc (2,j,i)
      end do
    end do

    write(nout) nvor_b

    do i = 1, nvor_b
      write(nout) vor_b    (1,i), vor_b    (2,i)
      write(nout) vor_vb (1,1,i), vor_vb (2,1,i)
      write(nout) vor_vb (1,2,i), vor_vb (2,2,i)
      write(nout) vor_vb (1,3,i), vor_vb (2,3,i)
      write(nout) cir_b      (i)
      write(nout) cor_b      (i)
      write(nout) blob_id    (i)
    end do

    write(nout) nvor_s

    do i = 1, nvor_s
      write(nout) vor_sc   (1,i), vor_sc   (2,i)
      write(nout) vor_s  (1,1,i), vor_s  (2,1,i)
      write(nout) vor_s  (1,2,i), vor_s  (2,2,i)
      write(nout) vor_vs (1,1,i), vor_vs (2,1,i)
      write(nout) vor_vs (1,2,i), vor_vs (2,2,i)
      write(nout) vor_vs (1,3,i), vor_vs (2,3,i)
      write(nout) cir_s      (i)
      write(nout) cor_s      (i)
      write(nout) vor_sr     (i)
      write(nout) sheet_id   (i)
    end do

   close(nout)

   call cp_out

  return
  end


! ************************************************************************** 
     subroutine cp_out
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

   character*40 nfile
   character*40 fout_h
   character*40 fout_t

   fout_h = '../data_out/p_'
   fout_t = '_step.dat'
   write(nfile,'(a14,i5.5,a9)') fout_h, nt, fout_t

   nout = 16
   open(nout,file=nfile,status='unknown',form='unformatted')

    write(nout) npanel

    do i = 1, npanel
      write(nout) poi_c(1,i),poi_c(2,i),cp_t(i),cf(i)
    end do

   close(nout)

  return
  end


