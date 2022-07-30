! ************************************************************************** 
!  convert from s_ data to vtk format for paraview                           
!   (vortex elements distribution at each time step can be displayed)        
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS CODE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING       
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 

  include './includes/modules.f90'
  include './includes/d_biot_p.f90'
  include './includes/c_solve_matrix_pressure.f90'
  include './LAPACK_library/dgemm.f'
  include './LAPACK_library/dger.f'
  include './LAPACK_library/dgesv.f'
  include './LAPACK_library/dgetf2.f'
  include './LAPACK_library/dgetrf.f'
  include './LAPACK_library/dgetrs.f'
  include './LAPACK_library/dlamch.f'
  include './LAPACK_library/dlaswp.f'
  include './LAPACK_library/dscal.f'
  include './LAPACK_library/dswap.f'
  include './LAPACK_library/dtrsm.f'
  include './LAPACK_library/idamax.f'
  include './LAPACK_library/ieeeck.f'
  include './LAPACK_library/ilaenv.f'
  include './LAPACK_library/iparmq.f'
  include './LAPACK_library/lsame.f'
  include './LAPACK_library/xerbla.f'

 use global
 use parameter
 use math_constant
 implicit none

 character(40) :: nfile
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

   nfile = '../data_out/velo-modified-cd-cl.dat'

    open  (11,file=nfile,status='unknown',form='formatted')
    write(11, *) 'ID,  nt,     time,      velocity,      modified_cx,     modified_cy'

   nfile = '../data_out/modified-cd-cl.dat'

    open  (12,file=nfile,status='unknown',form='formatted')
    write(12, *) 'ID,  nt,     time,      velocity,      modified_cx,     modified_cy'

!====   modify by k.k 2014/11/13   ===
    nfile = '../data_out/k_modified-cd-cl.dat'
    open  (13,file=nfile,status='unknown',form='formatted')
    write(13, *) 'ID,  nt,     time,      velocity,      modified_cpx,     modified_cfx      modified_cpy,     modified_cfy'
!=====================================

   do i = n_start, n_end, n_int
    call read_in  (i)
    call data_out (i)
   end do

   close(11)
   close(12)

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
   integer :: i

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
    read(10) dh_b
    read(10) dhdt
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
    read(10) edge_1  (i)
    read(10) edge_2  (i)
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
    read(10) panel_n_1(i)
    read(10) panel_n_2(i)
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
         cir_b(nvor_b) = cir_s     (i)
         cor_b(nvor_b) = dsqrt(vor_sr(i)*2.0*cor_s(i)/pi)
       blob_id(nvor_b) = sheet_id  (i)
   end do

   call solve_matrix_pressure


  return
  end

! ************************************************************************** 
     subroutine data_out (int)
! ************************************************************************** 

  use global
  use parameter
  use math_constant
  implicit none

   integer :: n2
   integer :: n4
   integer :: nc
   integer :: n_check
   integer :: i, j, k
   integer :: int
   real(8) :: vtt
   real(8) :: cx, cy, cz
   real(8) :: angle
   real(8) :: u2, v2
   real(8) :: ux, uy
   real(8) :: x1, y1
   real(8) :: ct, coef_cp

   real(8), dimension(:,:,:), allocatable :: cc

   allocate( cc(n_id, 2, 3) )

!===========================================================

   do i = 1, npanel
    cf(i) = 0.0d0
   end do

   do i = 1, npanel
    cf(i) = (uw(1,i)*poi_n(2,i) - uw(2,i)*poi_n(1,i)) /dh / re * 2.0d0
   end do


   angle = -theta/180.0d0*pi

   cx = 0.0d0
   cy = 0.0d0
   ct = 0.0d0

   do k = 1, nw
   do i = 1, 2
   do j = 1, 3
    cc(k,i,j) = 0.0d0
   end do
   end do
   end do

   coef_cp = 2.0d0

   do k = 1, nw

   do i = 1, npanel

!=================================
    if (panel_id(i) .eq. k) then  
!=================================

!    ////////// à≥óÕåWêîÇÃåvéZ //////////*

     x1 = poi_c(1,i)
     y1 = poi_c(2,i)

     ux = uw(1,i)
     uy = uw(2,i)
     u2 = ux*ux + uy*uy

     uinf = uniform_velocity
     vinf = 0.0d0

     v2 = ( uinf - vm(1,1,i) )**2 + ( vinf - vm(2,1,i) )**2

     cp_t(i) = coef_cp * ( ph(i) - 0.5d0*u2 + 0.5d0*v2 )

     rot1(i) =  poi_n(1,i)*dcos(angle) - poi_n(2,i)*dsin(angle)
     rot2(i) =  poi_n(1,i)*dsin(angle) + poi_n(2,i)*dcos(angle)

     rot3(i) =  poi_c(1,i)*poi_n(2,i)  - poi_c(2,i)*poi_n(1,i)
     rot4(i) = -poi_c(1,i)*poi_n(1,i)  + poi_c(2,i)*poi_n(2,i)

!    ////////// ó¨ëÃóÕ (à≥óÕê¨ï™) //////////*

     cp_x(i) = -cp_t(i)*rot1(i)*ds(i)
     cp_y(i) = -cp_t(i)*rot2(i)*ds(i)

     cc(k,1,1) = cc(k,1,1) - cp_t(i)*rot1(i)*ds(i)
     cc(k,1,2) = cc(k,1,2) - cp_t(i)*rot2(i)*ds(i)
     cc(k,1,3) = cc(k,1,3) + cp_t(i)*rot3(i)*ds(i)

!    ////////// ó¨ëÃóÕ (ñÄéCê¨ï™) //////////*

     cf_x(i) =  cf(i)*rot2(i)*ds(i)
     cf_y(i) = -cf(i)*rot1(i)*ds(i)

     cc(k,2,1) = cc(k,2,1) + cf(i)*rot2(i)*ds(i)
     cc(k,2,2) = cc(k,2,2) - cf(i)*rot1(i)*ds(i)
     cc(k,2,3) = cc(k,2,3) - cf(i)*rot4(i)*ds(i)

!=================================
    end if
!=================================

    end do

     cx = cx + cc(k,1,1) + cc(k,2,1)
     cy = cy + cc(k,1,2) + cc(k,2,2)
     ct = ct + cc(k,1,3) + cc(k,2,3)

   end do


!===========================================================

   do k = 1, nw
    write(11,'(i3,i6,f12.7,3e16.8)') k, nt, time, uniform_velocity-velocity,  &
                         (cc(k,1,1)+cc(k,2,1))/(uniform_velocity-velocity)  &
                                              /(uniform_velocity-velocity), &
                         (cc(k,1,2)+cc(k,2,2))/(uniform_velocity-velocity)  &
                                              /(uniform_velocity-velocity)
   end do


   do k = 1, nw
    write(12,'(i3,i6,f12.7,3e16.8)') k, nt, time, uniform_velocity-velocity,  &
                                     (cc(k,1,1)+cc(k,2,1)), &
                                     (cc(k,1,2)+cc(k,2,2))
   end do

  do k = 1, nw
    write(13,'(i3,i6,f12.7,5e16.8)') k, nt, time, uniform_velocity-velocity,  &
                                      cc(k,1,1)/(uniform_velocity-velocity)   &
                                               /(uniform_velocity-velocity),  &
                                      cc(k,2,1)/(uniform_velocity-velocity)   &
                                               /(uniform_velocity-velocity),  & ! cx
                                      cc(k,1,2)/(uniform_velocity-velocity)   &
                                               /(uniform_velocity-velocity),  &
                                      cc(k,2,2)/(uniform_velocity-velocity)   &
                                               /(uniform_velocity-velocity)    ! cy
   end do

  deallocate( cc )

  return
  end

