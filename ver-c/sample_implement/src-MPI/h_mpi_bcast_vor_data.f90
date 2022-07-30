! ************************************************************************** 
     subroutine mpi_bcast_vor_data
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

   nb = nvor_b
   ns = nvor_s

   call MPI_BCAST( blob_id   (1),  nb,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( cir_b     (1),  nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( cor_b     (1),  nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_b   (1,1),2*nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_vb(1,1,1),6*nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )

   call MPI_BCAST( sheet_id  (1),  ns,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( cir_s     (1),  ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( cor_s     (1),  ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_sr    (1),  ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_sc  (1,1),2*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_s (1,1,1),4*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )
   call MPI_BCAST( vor_vs(1,1,1),6*ns,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr )

  return
  end


