! ************************************************************************** 
      subroutine force
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! CALCULATION OF FLUID FORCE COFFICIENT                                      
!   cp_t : pressure cofficient                                               
!   cx : drag cofficient, cy : lift cofficient, ct : moment cofficient       
!============================================================================

  use global
  use parameter
  use math_constant
  implicit none
  include 'inc_2d'

   character(50) :: nfile1
   character(50) :: nfile2
   character(50) :: nfile3
   character(50) :: fin_h
   character(50) :: fin_m
   character(50) :: fin_t

   double precision, allocatable :: cc(:,:,:)

   allocate( cc(n_id,2,3) )

   call fric

   angle = -theta/180.0d0*pi

   cx = 0.0d0
   cy = 0.0d0
   ct = 0.0d0

   cc(1:nw,1:2,1:3) = 0.0d0

!===============
   do k = 1, nw 
!===============

   do i = 1, npanel

! ==================================
     if ( panel_id(i) .eq. k ) then 
! ==================================

     x1 = poi_c(1,i)
     y1 = poi_c(2,i)

     ux = lay_vc(1,nlayer(i),i)
     uy = lay_vc(2,nlayer(i),i)
     u2 = ux**2 + uy**2

     v2 = ( uinf - vm(1,1,i) )**2 + ( vinf - vm(2,1,i) )**2

     cp_t(i) = 2.0d0 * ( ph(i) - 0.5d0*u2 + 0.5d0*v2 )

     rot1(i) =  poi_n(1,i)*dcos(angle) - poi_n(2,i)*dsin(angle)
     rot2(i) =  poi_n(1,i)*dsin(angle) + poi_n(2,i)*dcos(angle)
     rot3(i) =  poi_c(1,i)*poi_n(2,i)  - poi_c(2,i)*poi_n(1,i)
     rot4(i) = -poi_c(1,i)*poi_n(1,i)  + poi_c(2,i)*poi_n(2,i)

! ////////// Cal. pressure cofficient ///////////

     cp_x(i) = -cp_t(i)*rot1(i)*ds(i)
     cp_y(i) = -cp_t(i)*rot2(i)*ds(i)

     cc(k,1,1) = cc(k,1,1) - cp_t(i)*rot1(i)*ds(i)
     cc(k,1,2) = cc(k,1,2) - cp_t(i)*rot2(i)*ds(i)
     cc(k,1,3) = cc(k,1,3) + cp_t(i)*rot3(i)*ds(i)

! ////////// Cal. friction cofficient ///////////

     cf_x(i) = -cf(i)*rot2(i)*ds(i)	! hattori
     cf_y(i) =  cf(i)*rot1(i)*ds(i)	! hattori

     cc(k,2,1) = cc(k,2,1) - cf(i)*rot2(i)*ds(i)	! hattori
     cc(k,2,2) = cc(k,2,2) + cf(i)*rot1(i)*ds(i)	! hattori
     cc(k,2,3) = cc(k,2,3) + cf(i)*rot4(i)*ds(i)	! hattori

! ==========
     end if 
! ==========

   end do

   cx = cx + cc(k,1,1) + cc(k,2,1)
   cy = cy + cc(k,1,2) + cc(k,2,2)
   ct = ct + cc(k,1,3) + cc(k,2,3)

!=========
   end do 
!=========

   write(*, *) ' '
   write(*, *) '  ID,       cx,       cy,       ct   (pressure force)'

   do k = 1, nw
     write(*,'(i5,1x,3f10.6)') k, cc(k,1,1), cc(k,1,2), cc(k,1,3)
   end do

   write(*, *) ' '
   write(*, *) '  ID,       cx,       cy,       ct   (friction force)'

   do k = 1, nw
     write(*,'(i5,1x,3f10.6)') k, cc(k,2,1), cc(k,2,2), cc(k,2,3)
   end do

   write(*, *) ' '
   write(*, *) '  ID,       cx,       cy,       ct   (total force)'

   do k = 1, nw
     write(*,'(i5,1x,3f10.6)') k, cc(k,1,1)+cc(k,2,1), cc(k,1,2)+cc(k,2,2), cc(k,1,3)+cc(k,2,3)
   end do

!===============
   do k = 1, nw 
!===============

   write(45,'(e16.8,i5,e16.8,1x,e16.8,1x,e16.8,1x,e16.4,1x,e16.4,1x,e16.4)') &
     time, k, cc(k,1,1), cc(k,1,2), cc(k,1,3), &
           parts_position(1,k), &
           parts_position(2,k), &
           parts_position(3,k)
   write(46,'(e16.8,i5,e16.8,1x,e16.8,1x,e16.8,1x,e16.4,1x,e16.4,1x,e16.4)') &
     time, k, cc(k,2,1), cc(k,2,2), cc(k,2,3), &
           parts_position(1,k), &
           parts_position(2,k), &
           parts_position(3,k)
   write(47,'(e16.8,i5,e16.8,1x,e16.8,1x,e16.8,1x,e16.4,1x,e16.4,1x,e16.4)') &
     time, k, cc(k,1,1)+cc(k,2,1), cc(k,1,2)+cc(k,2,2), cc(k,1,3)+cc(k,2,3), &
           parts_position(1,k), &
           parts_position(2,k), &
           parts_position(3,k)

!=========
   end do 
!=========

   write(48,'(e16.8,1x,e16.8,1x,e16.8,1x,e16.8)') time, cx, cy, ct

   fin_h = '../data_out/panel_force_'
   fin_t = '_step.dat'

   write(nfile1,'(a24,i5.5,a9)') fin_h,nt,fin_t
   open(49,file=nfile1,status='unknown',form='formatted')

! ===============
    do k = 1, nw 
! ===============

    write(49,*) 'Panel_ID=',k

    do i = 1, npanel

      if ( panel_id(i) .eq. k ) then 
        cp_f =  cp_x(i)*poi_n(1,i) + cp_y(i)*poi_n(2,i)
        cf_f =  cf_x(i)*( poi(1,2,i) - poi(1,1,i) ) + cf_y(i)*( poi(2,2,i) - poi(2,1,i) )
        write(49,'(e16.8,1x,e16.8)') cp_f, cf_f
      end if 

    end do

    write(49,*) ''

! =========
    end do 
! =========

   close (49)

   fin_h = '../data_out/panel_position_'
   fin_t = '_step.dat'

   write(nfile3,'(a27,i5.5,a9)') fin_h,nt,fin_t
   open(51,file=nfile3,status='unknown',form='formatted')

! ===============
    do k = 1, nw 
! ===============
    write(51,'(i2,e16.8,1x,e16.8,1x,e16.8)') &
    k, parts_position(1,k), parts_position(2,k), parts_position(3,k)
! =========
    end do 
! =========

   close (51)

   deallocate( cc )

  return
  end


! ************************************************************************** 
     subroutine fric
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

   cf(1:npanel) = 0.0d0

   do i = 1, npanel

     ux = lay_vc(1,1,i)
     uy = lay_vc(2,1,i)
     omega_s = ( uy*poi_n(1,i) - ux*poi_n(2,i) ) / lay_h(i)
     cf(i)   = omega_s / re * 2.0d0

   end do

  return
  end


