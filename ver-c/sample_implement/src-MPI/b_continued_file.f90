! ************************************************************************** 
     subroutine continued_file
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

   nout = 10
   open(nout,file=file2,status='unknown',form='unformatted')

    read(nout) re, dt
    read(nout) theta
    read(nout) omz
    read(nout) nt, time
    read(nout) dh
    read(nout) vn_vis
    read(nout) ds_t
    read(nout) dt_base
    read(nout) accelerate
    read(nout) velocity
    read(nout) velocity_base
    read(nout) uniform_velocity        ! add by K.F. (2012/08/13)
    read(nout) uniform_velocity_base   ! add by K.F. (2012/08/13)
    read(nout) u_accelerate            ! add by K.F. (2012/08/13)

    read(nout) nw
    do i = 1, nw
      read(nout) parts_position (1,i)
      read(nout) parts_position (2,i)
      read(nout) parts_position (3,i)
      read(nout) pgw   (1,i)
      read(nout) pgw   (2,i)
      read(nout) pgw   (3,i)
      read(nout) phase   (i)
      read(nout) cir_r   (i)
      read(nout) edge  (1,i)
      read(nout) edge  (2,i)
    end do

    read(nout) npanel
    do i = 1, npanel
      read(nout) q        (i)
      read(nout) ph       (i)
      read(nout) ds       (i)
      read(nout) cir_g    (i)
      read(nout) uw     (1,i), uw     (2,i)
      read(nout) poi  (1,1,i), poi  (2,1,i)
      read(nout) poi  (1,2,i), poi  (2,2,i)
      read(nout) poi_b(1,1,i), poi_b(2,1,i)
      read(nout) poi_b(1,2,i), poi_b(2,2,i)
      read(nout) poi_c  (1,i), poi_c  (2,i)
      read(nout) poi_n  (1,i), poi_n  (2,i)
      read(nout) vm   (1,1,i), vm   (2,1,i)
      read(nout) vm   (1,2,i), vm   (2,2,i)
      read(nout) am   (1,1,i), am   (2,1,i)
      read(nout) am   (1,2,i), am   (2,2,i)
      read(nout) panel_id (i)
      read(nout) panel_n(1,i)
      read(nout) panel_n(2,i)
    end do

    do i = 1, npanel
      read(nout) nlayer   (i)
      read(nout) lay_h    (i)
      read(nout) ypls     (i)
      do j = 1, nlayer(i)
        read(nout) lay_x  (1,j,i), lay_x  (2,j,i)
        read(nout) lay_y  (1,j,i), lay_y  (2,j,i)
        read(nout) lay_c  (1,j,i), lay_c  (2,j,i)
        read(nout) lay_vx (1,j,i), lay_vx (2,j,i)
        read(nout) lay_vy (1,j,i), lay_vy (2,j,i)
        read(nout) lay_vc (1,j,i), lay_vc (2,j,i)
      end do
    end do

    read(nout) nvor_b
    do i = 1, nvor_b
      read(nout) vor_b    (1,i), vor_b    (2,i)
      read(nout) vor_vb (1,1,i), vor_vb (2,1,i)
      read(nout) vor_vb (1,2,i), vor_vb (2,2,i)
      read(nout) vor_vb (1,3,i), vor_vb (2,3,i)
      read(nout) cir_b      (i)
      read(nout) cor_b      (i)
      read(nout) blob_id    (i)
    end do

    read(nout) nvor_s
    do i = 1, nvor_s
      read(nout) vor_sc   (1,i), vor_sc   (2,i)
      read(nout) vor_s  (1,1,i), vor_s  (2,1,i)
      read(nout) vor_s  (1,2,i), vor_s  (2,2,i)
      read(nout) vor_vs (1,1,i), vor_vs (2,1,i)
      read(nout) vor_vs (1,2,i), vor_vs (2,2,i)
      read(nout) vor_vs (1,3,i), vor_vs (2,3,i)
      read(nout) cir_s      (i)
      read(nout) cor_s      (i)
      read(nout) vor_sr     (i)
      read(nout) sheet_id   (i)
    end do

   close(nout)

  return
  end


