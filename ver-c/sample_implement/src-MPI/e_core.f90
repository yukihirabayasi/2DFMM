! ************************************************************************** 
     subroutine core
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

   corc = 2.2418d0

   do i = 1, nvor_b
     cor_b(i) = dsqrt( cor_b(i)**2 + (corc**2)*dt/re )
   end do

   do i = 1, nvor_s
     cor_s(i) = dsqrt( cor_s(i)**2 + (corc**2)*dt/re )
   end do

  return
  end


