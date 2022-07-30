! ************************************************************************** 
     subroutine mpi_bcast_panel_update_data
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

   np = npanel

   call MPI_BCAST( cir_r(1),  nw,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( phase(1),  nw,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( pgw(1,1),3*nw,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )

   call MPI_BCAST( q        (1),  np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( ph       (1),  np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( poi  (1,1,1),4*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( poi_b(1,1,1),4*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( uw     (1,1),2*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( cir_g    (1),  np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( poi_c  (1,1),2*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( ds       (1),  np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( poi_n  (1,1),2*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vm   (1,1,1),4*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( am   (1,1,1),4*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )

  return
  end


