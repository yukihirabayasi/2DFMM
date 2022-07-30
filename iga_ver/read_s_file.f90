include './includes/modules.f90'

module m_read_s_file
  implicit none
contains
  subroutine read_s_file(vor)
    use m_QuadTree
  
    use global
    use parameter
    use math_constant
  
    type(Vortex),intent(out),allocatable :: vor(:)
  
    character(50), save :: f_inp
    character(50), save :: f_debug
    integer nf_inp
    integer nf_debug
    integer i,j

    print *,"- read s_file"

    nf_inp = 10 
    f_inp = "./s_00344_step.dat" 
    nf_debug = 300
    f_debug = "./vor.dat"
  
    open(nf_inp,file=f_inp,status='old',form='unformatted')
    read(nf_inp) re, dt
    read(nf_inp) theta
    read(nf_inp) omz
    read(nf_inp) nt, time
    read(nf_inp) dh
    read(nf_inp) ds_t
    read(nf_inp) i_mirror
    read(nf_inp) d_mirror
    read(nf_inp) pheta
    read(nf_inp) dt_base
    read(nf_inp) accelerate
    read(nf_inp) velocity
    read(nf_inp) velocity_base
    read(nf_inp) uniform_velocity        ! add by K.F. (2012/08/13)
    read(nf_inp) uniform_velocity_base   ! add by K.F. (2012/08/13)
    read(nf_inp) u_accelerate            ! add by K.F. (2012/08/13)
  
    read(nf_inp) nw
    do i = 1, nw
      read(nf_inp) parts_position (1,i)
      read(nf_inp) parts_position (2,i)
      read(nf_inp) parts_position (3,i)
      read(nf_inp) cir_r   (i)
      read(nf_inp) edge  (1,i)
      read(nf_inp) edge  (2,i)
    end do
  
    read(nf_inp) npanel
    do i = 1, npanel
      read(nf_inp) q         (i)
      read(nf_inp) ph        (i)
      read(nf_inp) ds        (i)
      read(nf_inp) cir_g     (i)
      read(nf_inp) uw      (1,i), uw      (2,i)
      read(nf_inp) poi   (1,1,i), poi   (2,1,i)
      read(nf_inp) poi   (1,2,i), poi   (2,2,i)
      read(nf_inp) poi_b (1,1,i), poi_b (2,1,i)
      read(nf_inp) poi_b (1,2,i), poi_b (2,2,i)
      read(nf_inp) poi_c   (1,i), poi_c   (2,i)
      read(nf_inp) poi_n   (1,i), poi_n   (2,i)
      read(nf_inp) vm    (1,1,i), vm    (2,1,i)
      read(nf_inp) vm    (1,2,i), vm    (2,2,i)
      read(nf_inp) am    (1,1,i), am    (2,1,i)
      read(nf_inp) am    (1,2,i), am    (2,2,i)
      read(nf_inp) panel_n (1,i), panel_n (2,i)
      read(nf_inp) panel_id  (i)
    end do
  
    do i = 1, npanel
      read(nf_inp) nlayer         (i)
      read(nf_inp) lay_height     (i)
      read(nf_inp) y_plus         (i)
      do j = 1, nlayer(i)
      read(nf_inp) lay_x      (1,j,i), lay_x      (2,j,i)
      read(nf_inp) lay_y      (1,j,i), lay_y      (2,j,i)
      read(nf_inp) lay_center (1,j,i), lay_center (2,j,i)
      read(nf_inp) lay_velo   (1,j,i), lay_velo   (2,j,i)
      end do
    end do
  
    read(nf_inp) nvor_b
    allocate(vor(nvor_b))
    do i = 1, nvor_b
      read(nf_inp) vor_b    (1,i), vor_b    (2,i)
      read(nf_inp) vor_vb (1,1,i), vor_vb (2,1,i)
      read(nf_inp) vor_vb (1,2,i), vor_vb (2,2,i)
      read(nf_inp) vor_vb (1,3,i), vor_vb (2,3,i)
      read(nf_inp) cir_b      (i)
      read(nf_inp) cor_b      (i)
      read(nf_inp) blob_id    (i)
      vor(i)%pos(1) = vor_b(1,i)
      vor(i)%pos(2) = vor_b(2,i)
      vor(i)%cir = cir_b(i)
      vor(i)%cor = cor_b(i)
    end do

    read(nf_inp) nvor_s
    do i = 1, nvor_s
      read(nf_inp) vor_sc   (1,i), vor_sc   (2,i)
      read(nf_inp) vor_s  (1,1,i), vor_s  (2,1,i)
      read(nf_inp) vor_s  (1,2,i), vor_s  (2,2,i)
      read(nf_inp) vor_vs (1,1,i), vor_vs (2,1,i)
      read(nf_inp) vor_vs (1,2,i), vor_vs (2,2,i)
      read(nf_inp) vor_vs (1,3,i), vor_vs (2,3,i)
      read(nf_inp) cir_s      (i)
      read(nf_inp) cor_s      (i)
      read(nf_inp) vor_sr     (i)
      read(nf_inp) sheet_id   (i)
    end do
  
    close(nf_inp)
    
    open(nf_debug,file=f_debug,status="replace")
    do i = 1, nvor_b
      write(nf_debug,*) vor(i)%pos(1),",",vor(i)%pos(2)
    end do
    close(nf_debug)
  
  end subroutine read_s_file
end module
