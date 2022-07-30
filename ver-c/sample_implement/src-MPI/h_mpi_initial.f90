! ************************************************************************** 
     subroutine mpi_initial
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

   call MPI_INIT     ( ierr )
   call MPI_COMM_RANK( MPI_COMM_WORLD,my_rank,  ierr )
   call MPI_COMM_SIZE( MPI_COMM_WORLD,num_procs,ierr )

  return
  end


