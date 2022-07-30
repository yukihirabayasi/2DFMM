! ************************************************************************** 
      subroutine move
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

   do i = 1, npanel
    poi(1,1,i) = poi_b(1,1,i)
    poi(2,1,i) = poi_b(2,1,i)
    poi(1,2,i) = poi_b(1,2,i)
    poi(2,2,i) = poi_b(2,2,i)
    vm (1,2,i) = vm   (1,1,i)
    vm (2,2,i) = vm   (2,1,i)
    vm (1,1,i) = 0.0d0
    vm (2,1,i) = 0.0d0
    am (1,2,i) = am   (1,1,i)
    am (2,2,i) = am   (2,1,i)
    am (1,1,i) = 0.0d0
    am (2,1,i) = 0.0d0
   end do

   return
   end


