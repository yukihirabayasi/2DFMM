! ************************************************************************** 
     subroutine panel
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

!////////// Panel Center //////////*

   do i = 1, npanel
     poi_c(1,i) = 0.5d0*( poi(1,1,i) + poi(1,2,i) )
     poi_c(2,i) = 0.5d0*( poi(2,1,i) + poi(2,2,i) )
   end do

!////////// panel normal vector and area //////////*

   do i = 1, npanel
     rx = poi(1,2,i) - poi(1,1,i)
     ry = poi(2,2,i) - poi(2,1,i)
     ds(i) = dsqrt(rx**2 + ry**2)
     poi_n(1,i) =  ry / ds(i)
     poi_n(2,i) = -rx / ds(i)
   end do

!///////////// panel area /////////////*

   ds_t = 0.0d0

   do i = 1, npanel
     ds_t = ds_t + ds(i)
   end do

! ------- add by kuji (2015/03/19) -----

!//////////// edge_cir ///////////////*

   do j = 1, nw
     edge_cir(1,j) = 0.0d0
     edge_cir(2,j) = 0.0d0
   end do

! -------------------------------------
  return
  end


! ************************************************************************** 
     subroutine panel_height
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

   if ( i_panel_height_check .eq. 1 ) then
     dh = c_h * 5.0d0 * dsqrt( 1.0d0/re )
   else
     dh = c_h * ds_t / npanel
   end if

   write(*,*)              ''
   write(*,'(2x,A,f10.6)') 'dh    =',dh
   write(*,*)              ''

  return
  end


