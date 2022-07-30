include './tree/m_tree_bindc.f90'

! ************************************************************************** 
      subroutine get_velo_mpi
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
! Calculation of convective velocity of vortex elements.                     
!============================================================================

  use global
  use parameter
  use m_tree_bindc
  use math_constant
  
  implicit none
  include 'inc_2d'
  include 'mpif.h'

   integer istatus (MPI_STATUS_SIZE)
   integer :: mpi_start
   integer :: mpi_end

!**************************************
!---// velocity for blob elements //---
!**************************************

   call mpi_load_balance( nvor_b,my_rank,mpi_start,mpi_end )
   call make_tree(nvor_b,vor_b,cir_b,cor_b)

!   write (*,'(A,i3,A,i8,A,i8,A,i6)') '   blob  (rank ',my_rank,')' &
!  & ,mpi_start,'  to ',mpi_end,'   ==>   allocated : ',mpi_end - mpi_start + 1

   do i = mpi_start, mpi_end
     x1 = vor_b(1,i)
     y1 = vor_b(2,i)
     call indus_s( x1,y1,u1,v1 )
     call indus_v( x1,y1,u2,v2 )
     call biot_tree ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )
     vmx = uinf + ( -omz * y1 )
     vmy = vinf + (  omz * x1 )
     vor_vb(1,1,i) = vmx + ( u1 + u2 + u3 + u4 )
     vor_vb(2,1,i) = vmy + ( v1 + v2 + v3 + v4 )
   end do

!***************************************
!---// velocity for sheet elements //---
!***************************************

   call mpi_load_balance( nvor_s,my_rank,mpi_start,mpi_end )

!   write (*,'(A,i3,A,i8,A,i8,A,i6)') '   sheet (rank ',my_rank,')' &
!  & ,mpi_start,'  to ',mpi_end,'   ==>   allocated : ',mpi_end - mpi_start + 1

   do i = mpi_start, mpi_end
     x1 = vor_sc(1,i)
     y1 = vor_sc(2,i)
     call indus_s( x1,y1,u1,v1 )
     call indus_v( x1,y1,u2,v2 )
     call biot_tree ( x1,y1,u3,v3 )
     call biot_s ( x1,y1,u4,v4 )
     vmx = uinf + ( -omz * y1 )
     vmy = vinf + (  omz * x1 )
     vor_vs(1,1,i) = vmx + ( u1 + u2 + u3 + u4 )
     vor_vs(2,1,i) = vmy + ( v1 + v2 + v3 + v4 )
   end do

!****************************
!---// send and receive //---
!****************************

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

!*********************
!--- blob elements ---
!*********************
! ============================
   if ( my_rank .ne. 0 ) then 
! ============================
   call mpi_load_balance ( nvor_b,my_rank,  mpi_start,mpi_end )
   call mpi_send_vor_velo( vor_vb,nvor_b, 0,mpi_start,mpi_end )
! =================================
   else if ( my_rank .eq. 0 ) then 
! =================================
   do i = 1, num_procs-1
     call mpi_load_balance ( nvor_b,       i,mpi_start,mpi_end )
     call mpi_recv_vor_velo( vor_vb,nvor_b,i,mpi_start,mpi_end )
   end do
! ========
   end if 
! ========

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

!**********************
!--- sheet elements ---
!**********************
! ============================
   if ( my_rank .ne. 0 ) then 
! ============================
   call mpi_load_balance  ( nvor_s,my_rank,  mpi_start,mpi_end )
   call mpi_send_vor_velo ( vor_vs,nvor_s, 0,mpi_start,mpi_end )
! =================================
   else if ( my_rank .eq. 0 ) then 
! =================================
   do i = 1, num_procs-1
     call mpi_load_balance ( nvor_s,       i,mpi_start,mpi_end )
     call mpi_recv_vor_velo( vor_vs,nvor_s,i,mpi_start,mpi_end )
   end do
! ========
   end if 
! ========

   call MPI_BARRIER( MPI_COMM_WORLD,ierr )

  return
  end


