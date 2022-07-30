! ************************************************************************** 
     subroutine mpi_load_balance(n_vor_mpi,n_rank,n_start,n_end)
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
  include 'mpif.h'

   integer :: n_vor_mpi
   integer :: n_div
   integer :: n_sub
   integer :: n_rank
   integer :: n_start
   integer :: n_end

   n_sub = mod(n_vor_mpi,num_procs)
   n_div = int(n_vor_mpi/num_procs)

   do i = 0, n_rank
     if ( i .eq. 0 ) then 
       n_start = 1
       if ( n_sub .ne. 0 ) n_end = n_div + 1
       if ( n_sub .eq. 0 ) n_end = n_div
     else
       n_start = n_end + 1
       if ( i .le. n_sub-1 ) then
         n_end = n_start + n_div
       else
         n_end = n_start + n_div - 1
       end if
     end if
   end do

  return
  end


