! ************************************************************************** 
     subroutine body
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

   double precision, allocatable :: xp(:)
   double precision, allocatable :: yp(:)

   allocate( xp(n_panel) )
   allocate( yp(n_panel) )

   open (12, file=file1, status='old', form='formatted')

    read (12,*) nw

    write(*,*)           ''
    write(*,*)           ''
    write(*,*)           '*** Panel data ***'
    write(*,*)           ' Panel data   : ',file1
    write(*,'(2x,A,i2)') 'Parts number : ', nw

    n0     = 0
    npanel = 0

! ===============
    do i = 1, nw 
! ===============

    read(12,*) n

    write(*,'(2x,A,i2,1x,A,i5,A)') '   ID : ',i,'  includes ',n, '   elements'
    do j = 1, n + 1
      read(12,*) xp(j), yp(j)
    end do

    do j = 1, n
      npanel = npanel + 1

      poi_b(1, 1, npanel) = xp(j)
      poi_b(2, 1, npanel) = yp(j)
      poi_b(1, 2, npanel) = xp(j+1)
      poi_b(2, 2, npanel) = yp(j+1)
      panel_id   (npanel) = i
      panel_n  (1,npanel) = npanel - 1
      panel_n  (2,npanel) = npanel + 1

      ! modified by KF on the 6th of August, 2015

      if ( j .eq. 1 .and. i_edge .eq. 0 ) edge(1,i) = npanel
      if ( j .eq. n .and. i_edge .eq. 0 ) edge(2,i) = npanel
      if ( j .eq. 1 .and. i_edge .eq. 1 ) edge(1,i) = 1
      if ( j .eq. n .and. i_edge .eq. 1 ) edge(2,i) = 1

      if ( j .eq. 1 .and. i_edge .eq. 0 ) panel_n(1,npanel) = n
      if ( j .eq. n .and. i_edge .eq. 0 ) panel_n(2,npanel) = 1
      if ( j .eq. 1 .and. i_edge .eq. 1 ) panel_n(1,npanel) = 1
      if ( j .eq. n .and. i_edge .eq. 1 ) panel_n(2,npanel) = n

      if ( nt .eq. 0 ) then
        cir_g   (npanel) = 0.0d0
        vm  (1:2,1,npanel) = 0.0d0
        am  (1:2,1,npanel) = 0.0d0
      end if

    end do

    if ( nt .eq. 0 ) then
      parts_position (1:3,i) = 0.0d0
      pgw            (1:3,i) = 0.0d0
      phase              (i) = 0.0d0
      cir_r              (i) = 0.0d0
    end if

    n0 = npanel

! =========
    end do 
! =========

   close (12)

   deallocate( xp )
   deallocate( yp )

  return
  end


