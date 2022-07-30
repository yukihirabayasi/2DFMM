! ************************************************************************** 
     subroutine mpi_send_vor_velo(bcell,n_vor,n_rank,n_start,n_end)          
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

   integer :: n_trans
   integer :: n_rank
   integer :: n_start
   integer :: n_end

   double precision :: bcell (2,3,n_vor)
   double precision :: m1    (n_sp)
   double precision :: m2    (n_sp)

   n_trans = int((n_end - n_start)/n_sp) + 1

!====================
   do j = 1, n_trans 
!====================

   i = n_start + n_sp*(j-1)

   if ( j .eq. n_trans ) then
     np = n_end - i + 1
   else
     np = n_sp
   end if

   do n = 1, np
     m1(n) = bcell(1,1,i)
     m2(n) = bcell(2,1,i)
     i = i + 1
   end do

   n1 = my_rank*100 + j*10 + nt
   n2 = my_rank*200 + j*20 + nt

   call MPI_SEND( m1(1),np,MPI_DOUBLE_PRECISION,n_rank,n1,MPI_COMM_WORLD,ierr )
   call MPI_SEND( m2(1),np,MPI_DOUBLE_PRECISION,n_rank,n2,MPI_COMM_WORLD,ierr )

!=========
   end do 
!=========

  return
  end


