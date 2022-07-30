! ************************************************************************** 
     subroutine velo_update
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

   if ( my_rank .eq. 0 .and. nt .eq. 0 ) then
     do i = 1, nvor_b
       vor_vb(1,2,i) = vor_vb(1,1,i)
       vor_vb(2,2,i) = vor_vb(2,1,i)
     end do
     do i = 1, nvor_s
       vor_vs(1,2,i) = vor_vs(1,1,i)
       vor_vs(2,2,i) = vor_vs(2,1,i)
     end do
   end if

   do i = 1, nvor_b
     vor_vb(1,3,i) = vor_vb(1,2,i)
     vor_vb(2,3,i) = vor_vb(2,2,i)
     vor_vb(1,2,i) = vor_vb(1,1,i)
     vor_vb(2,2,i) = vor_vb(2,1,i)
   end do

   do i = 1, nvor_s
     vor_vs(1,3,i) = vor_vs(1,2,i)
     vor_vs(2,3,i) = vor_vs(2,2,i)
     vor_vs(1,2,i) = vor_vs(1,1,i)
     vor_vs(2,2,i) = vor_vs(2,1,i)
   end do

  return
  end


