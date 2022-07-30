! ************************************************************************** 
      subroutine biot_tr ( xi,yi,ax,ay )
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

   ax = 0.0d0
   ay = 0.0d0

   i = root

 10 continue

   if ( i .ne. 0 ) then

   rx = xi - pos(1,i)
   ry = yi - pos(2,i)
   r  = dsqrt( rx**2 + ry**2 )

!============================
   if ( i .lt. n_cell ) then 
!============================

   if ( r .gt. 1.0d-6 ) then

     rv2 = 1.0d0 / r**2
     xai = r / scl(i)
     rv  = vec(i) * ( 1.0d0 - dexp( -xai**2 ) )

     ax = ax + ry*rv*rv2
     ay = ay - rx*rv*rv2

   end if

   i = next(i)

!=======
   else 
!=======

   if ( r**2 .ge. r_critical(i) ) then

     rv2 = 1.0d0 / r**2
     xai = r / scl(i)
     rv  = vec(i) * ( 1.0d0 - dexp( -xai**2 ) )

     ax = ax + ry*rv*rv2
     ay = ay - rx*rv*rv2

     i = next(i)

   else

     i = more(i)

   end if

!=========
   end if 
!=========

   goto 10

   end if

   return

   end


! ************************************************************************** 
      subroutine make_tree
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

   call send_vortex_tree

   if ( n_body .gt. 0 ) then
     call expand_box
     call load_tree
     call cell_mass
     call thread
   end if

   return

   end


! ************************************************************************** 
      subroutine send_vortex_tree
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

   n_body = 0

   do i = 1, nvor_b

! ************ normal ******************

     n_body = n_body + 1

     pos (1,n_body) =  vor_b(1,i)
     pos (2,n_body) =  vor_b(2,i)
     vel (1,n_body) =  vor_vb(1,1,i)
     vel (2,n_body) =  vor_vb(2,1,i)
     mas   (n_body) =  dabs( cir_b(i) ) * pitwo
     vec   (n_body) =  cir_b(i) * pitwo
     scl   (n_body) =  cor_b(i)

! ************ mirror (y-direction) ******************

     if ( i_mirror .eq. 1 ) then

     n_body = n_body + 1

     pos (1,n_body) =  vor_b(1,i)
     pos (2,n_body) =  2.0d0*d_mirror - vor_b(2,i)
     vel (1,n_body) =  vor_vb(1,1,i)
     vel (2,n_body) = -vor_vb(2,1,i)
     mas   (n_body) =  dabs( cir_b(i) ) * pitwo
     vec   (n_body) = -cir_b(i) * pitwo
     scl   (n_body) =  cor_b(i)

     end if

! ************ end mirror image ******************

   end do

   return

   end


! ************************************************************************** 
      subroutine expand_box
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Expand root area to hold all particles.                                   
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

  double precision :: xymax

   xymax = 0.0d0
   size  = 4.0d0

   do i = 1, 2
   do j = 1, n_body
     xymax = dmax1( xymax,dabs( pos(i,j) ) )
   end do
   end do

 30 continue

   if ( xymax .ge. size/2.0d0 ) then
     size = 2.0d0 * size
     goto 30
   end if

   return

   end


! ************************************************************************** 
      subroutine load_tree
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Load particles into the tree.                                             
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

  integer :: makecell

   ncell = 0

   root = makecell( )

   midpoint(1:2,root) = 0.0d0

   cellsize(root) = size

   do i = 1, n_body
     call load_body( i )
   end do

   return

   end


! ************************************************************************** 
      subroutine load_body ( i )
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

  integer :: c
  integer :: i0
  integer :: jindex
  integer :: makecell
  integer :: subindex

!  ////////// start j, jindex pair in correct subcell of root ///////////

   j = root

   jindex = subindex(i,j)

 20 continue

   if ( subcell(jindex,j) .ne. 0      ) then

   if ( subcell(jindex,j) .lt. n_cell ) then

     c = makecell( )

     do k = 1, 2
       if ( pos(k,i) .ge. midpoint(k,j) ) then
         midpoint(k,c) = midpoint(k,j) + cellsize(j) / 4.0d0
       else
         midpoint(k,c) = midpoint(k,j) - cellsize(j) / 4.0d0
       end if
     end do

     cellsize(c) = cellsize(j) / 2.0d0

     i0 = subcell(jindex,j)

     subcell(subindex(i0,c),c) = i0
     subcell(jindex,j) = c

   end if

   j  = subcell(jindex,j)
   jindex = subindex(i,j)

   goto 20

   end if

   subcell(jindex,j) = i

   return

   end


! ************************************************************************** 
      subroutine cell_mass
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Compute cell masses, c.m. positions, check tree structure,                
!  assign critical radius, and compute quadrupole moments.                   
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

  integer :: index(n_cell)

  double precision :: distance
  double precision :: position(2)

!  ////////// list cells in order of decreasing size ///////////

   call breadth_first_search( index )

!  ////////// loop processing cells from smallest to root ///////////

!======================
   do i = ncell, 1, -1 
!======================

   n1 = index(i)

   mas  (n1) = 0.0d0
   vec  (n1) = 0.0d0
   scl  (n1) = 0.0d0
   vel(1,n1) = 0.0d0
   vel(2,n1) = 0.0d0

   do j = 1, 2
     position(j) = 0.0d0
   end do

   do j = 1, 4

     n2 = subcell(j,n1)

     if ( n2 .ne. 0 ) then

       mas  (n1) = mas  (n1) + mas(n2)
       vec  (n1) = vec  (n1) + vec(n2)
       scl  (n1) = scl  (n1) + mas(n2) * scl  (n2)
       vel(1,n1) = vel(1,n1) + mas(n2) * vel(1,n2)
       vel(2,n1) = vel(2,n1) + mas(n2) * vel(2,n2)

       do k = 1, 2
         position(k) = position(k) + mas(n2) * pos(k,n2)
       end do

     end if

   end do

   scl  (n1) = scl  (n1) / mas(n1)
   vel(1,n1) = vel(1,n1) / mas(n1)
   vel(2,n1) = vel(2,n1) / mas(n1)

   do j = 1, 2
     position(j) = position(j) / mas(n1)
   end do

   distance = 0.0d0

   do j = 1, 2
     distance  = distance + ( position(j) - midpoint(j,n1) )**2
     pos(j,n1) = position(j)
   end do

   r_critical(n1) = ( cellsize(n1) / pheta + dsqrt( distance ) )**2

!=========
   end do 
!=========

   return

   end


! ************************************************************************** 
      subroutine breadth_first_search ( index )
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

  integer :: firstcell
  integer :: lastcell
  integer :: nextcell

  integer :: index(n_cell)

!  ////////// start scan with root as only active cell ///////////

   index(1) = root

   firstcell = 1
   lastcell  = 1

!  ////////// loop while active cells to process ///////////

 10 continue

   if ( firstcell .le. lastcell ) then

!  ////////// start counting active cells in next iteration ///////////

   nextcell = lastcell

!  ////////// loop over subcells of each active cell ///////////

   do i = 1, 4
   do j = firstcell, lastcell

!  ////////// add all cells on next level to active list ///////////

     if ( subcell(i,index(j)) .ge. n_cell ) then
       nextcell = nextcell + 1
       index(nextcell) = subcell(i,index(j))
     end if

   end do
   end do

!  ////////// advance first and last active cell indices ///////////

   firstcell = lastcell + 1
   lastcell  = nextcell

   goto 10

   end if

   return

   end


! ************************************************************************** 
      subroutine thread
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Check tree structure.                                                     
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

  integer :: sptr
  integer :: stack(0:1024)

!  ////////// push root cell onto stack ///////////

   sptr = 1

   stack(sptr) = root

!  ////////// loop while stuff on stack to process ///////////

 10 continue

   if ( sptr .gt. 0 ) then

!  ////////// pop node on top of stack ///////////

   i = stack(sptr)

   sptr = sptr - 1

!  ////////// set index of next node to visit ///////////

   next(i) = stack(sptr)

!  ////////// push descendents if i is a cell ///////////

   if ( i .ge. n_cell ) then

!  ////////// loop over existing descendents ///////////

     do j = 4, 1, -1

       if ( subcell(j,i) .ne. 0 ) then

!  ////////// push each descendent on stack ///////////

         sptr = sptr + 1
         stack(sptr) = subcell(j,i)

       end if

     end do

!  ////////// set index of first descendent ///////////

     more(i) = stack(sptr)

   end if

   goto 10

   end if

   return

   end


! ************************************************************************** 
      integer function makecell( )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Function to allocate a cell, returning its index.                         
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

   ncell    = ncell + 1
   makecell = ncell + n_vor_b

   subcell(1:4,makecell) = 0

   return

   end


! ************************************************************************** 
      integer function subindex ( i,j )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Compute subcell index for node i within cell j.                           
!============================================================================

 use global
 use parameter
 use math_constant
 implicit none
 include 'inc_2d'

   subindex = 1

   do k = 1, 2
     if ( pos(k,i) .ge. midpoint(k,j) ) subindex = subindex + 2**( 2-k )
   end do

   return

   end

