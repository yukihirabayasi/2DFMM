module m_morton_order
  interface
    subroutine getChild(p_index,c_index)
      integer,intent(in) :: p_index
      integer,intent(out) :: c_index(4)
    end subroutine getChild

    subroutine getParent(c_index,p_index)
      integer,intent(in) :: c_index
      integer,intent(out) :: p_index
    end subroutine getParent

    subroutine getLevel(index,level)
      integer,intent(in) :: index
      integer,intent(out) :: level
    end subroutine getLevel

    subroutine index2pos(globalIndex,pos)
      integer,intent(in) :: globalIndex
      integer,intent(out) :: pos(2)
    end subroutine index2pos

    subroutine pos2index(pos,level,index)
      integer,intent(in) :: pos(2)
      integer,intent(in) :: level
      integer,intent(out) :: index
    end subroutine pos2index

    subroutine testChild(p_index,c_index,p_level)
      integer,intent(in) :: p_index
      integer,intent(in) :: p_level
      integer,intent(out) :: c_index(4)
    end subroutine testChild

    subroutine get_neighbor_set(t_index,neighbor_set) 
      integer,intent(in) :: t_index
      integer,intent(out) :: neighbor_set(8)
    end subroutine get_neighbor_set
  end interface
end module m_morton_order
!contains
!-------------------------------------------------
!- subroutine
!-------------------------------------------------
subroutine getChild(p_index, c_index) 
  integer,intent(in) :: p_index
  integer,intent(out) :: c_index(4)
  integer :: level
  integer :: index
  integer :: i

  call getLevel(p_index,level)
  index = p_index - ((4**level+1)-1)/3 - 1
  index = ishft(index,2)
  do i = 0,3 
    c_index(i+1) = index+i
  end do
  c_index = c_index + 1 + (4**(level+1)-1)/3
end subroutine getChild

subroutine testChild(p_index, c_index, p_level) 
  integer,intent(in) :: p_index
  integer,intent(in) :: p_level
  integer,intent(out) :: c_index(4)
  integer :: level
  integer :: index
  integer :: i

!  call getLevel(p_index,level)
  index = p_index - ((4**p_level+1)-1)/3 - 1
  index = ishft(index,2)
  do i = 0,3 
    c_index(i+1) = index+i
  end do
  c_index = c_index + 1 + (4**(p_level+1)-1)/3
end subroutine testChild


subroutine getParent(c_index,p_index)
  integer,intent(in) :: c_index
  integer,intent(out) :: p_index
  integer :: c_level
  call getLevel(c_index,c_level) 
  p_index = ishft(c_index-(4**c_level-1)/3-1,-2)+(4**(c_level-1)-1)/3+1
end subroutine getParent


subroutine getLevel(index,level)
  integer,intent(in) :: index
  integer,intent(out) :: level ! nodes(index)%level
  level = aint(log(3.d0*real(index))/log(4.d0))
end subroutine getLevel



subroutine index2pos(globalIndex,pos)
! split index -> position x,y
! example) 1011 -> 01,11 -> 1,3
  integer,intent(in) :: globalIndex
  integer,intent(out) :: pos(2)
  integer :: nlevel = 15 ! index2pos overflows when it exceeds 16bits
  integer :: i, index_xy
  integer :: localIndex
  integer :: level

  call getLevel(globalIndex,level)
  localIndex = globalIndex - (4**level-1)/3 - 1
  do i = 0,2*nlevel
    index_xy = mod(i,2)+1 ! 0 or 1
    if (btest(localIndex,i)) then
      pos(index_xy) = ibset(pos(index_xy), i/2)
    else
      pos(index_xy) = ibclr(pos(index_xy), i/2)
    end if
  end do

end subroutine index2pos

subroutine pos2index(pos,level,index) 
! merge position x,y -> index
! example) 1,3 -> 01,11 -> 1011
  integer,intent(in) :: pos(2)
  integer,intent(in) :: level
  integer,intent(out) :: index
  integer :: i, index_xy
  integer :: nlevel = 15
  index = 0
  do i = 0,2*nlevel
    index_xy = mod(i,2)+1
    if (btest(pos(index_xy),i/2)) then
      index = ibset(index,i)
    else
      index = ibclr(index,i)
    end if
  end do

  index = index + (4**level-1)/3 + 1

end subroutine pos2index


subroutine get_neighbor_set(t_index,neighbor_set) 
! Assign -1 to the array when there is no target cell
  integer,intent(in) :: t_index
  integer,intent(out) :: neighbor_set(8)
  integer :: t_pos(2)
  integer :: dummy_pos(2)
  integer :: counter
  integer :: test
  integer :: index
  integer :: level,minIndex,maxIndex
  counter = 0

  !- Calc MaxIndex -
  call getLevel(t_index,level)
  minIndex = (4**level-1)/3 + 1
  maxIndex = (4**(level+1)-1)/3

  call index2pos(t_index,t_pos)

  do i = -1,1,1
    do j = -1,1,1

      if (i==0 .and. j==0) then
        cycle
      else
        counter = counter + 1
      end if

      if (t_pos(1)+i < 0 .or. t_pos(2)+j < 0) then
        neighbor_set(counter) = 0
        cycle
      end if

      dummy_pos(1) = t_pos(1)+i
      dummy_pos(2) = t_pos(2)+j
      call pos2index(dummy_pos,level,index)
      if (index < minIndex .or. maxIndex < index) then
        neighbor_set(counter) = 0
      else 
        neighbor_set(counter) = index
      end if
    end do
  end do
end subroutine get_neighbor_set


subroutine get_interaction_set(t_index,interaction_set)
  integer,intent(in) :: t_index
  integer,intent(out) :: interaction_set(27) 
  integer :: p_index,c_index(4)
  integer :: p_neighbor_set(8)
  integer :: i,j,k,counter
  integer :: exclusion_set(8)

  interaction_set(:) = 0
  counter = 0
  call getParent(t_index, p_index)
  call get_neighbor_set(p_index, p_neighbor_set)
  call get_neighbor_set(t_index, exclusion_set)
  do i = 1, 8
    if (p_neighbor_set(i) == 0) cycle
    call getChild(p_neighbor_set(i),c_index)
    do j = 1, 4
      do k = 1, 8
        if (c_index(j) == exclusion_set(k)) exit 
        if (k == 8) then
          counter = counter + 1
          interaction_set(counter) = c_index(j)
        end if
      end do
    end do
  end do

end subroutine get_interaction_set