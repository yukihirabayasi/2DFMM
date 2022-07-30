! ************************************************************************** 
     subroutine mpi_bcast_cal_update_data
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

   call MPI_BCAST( nt,    1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( nw,    1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( npanel,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( nvor_b,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( nvor_s,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( dt,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( time,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( uinf,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vinf,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )

  return
  end


