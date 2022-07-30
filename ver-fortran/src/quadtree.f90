module m_QuadTree
  implicit none
  type :: Node
    integer :: interaction_set(27)
    integer :: level
    integer :: id
    integer,allocatable :: vor_id(:)
    integer :: num_vor
    double precision :: center_pos(2)
    double precision :: region ! |x_min| = |x_max| = |y_min| = |y_max|
    complex(kind(0d0)),allocatable :: outer_coeffs(:)
    complex(kind(0d0)),allocatable :: inner_coeffs(:)
    logical :: is_leaf
  end type Node

  type :: QuadTree
    type(Node),allocatable :: nodes(:)! => null()
    integer :: nterms
    integer :: max_level
    double precision :: theta ! tolerance parameter
    integer :: level_hist(15) ! for debug
    integer :: t_level ! for biot_FMM
    double precision :: r_min ! for Tree
    integer :: stat ! for maxLevel loop
  end type QuadTree

  type :: Vortex
    double precision :: pos(2)
    double precision :: cir 
    double precision :: cor
  end type Vortex

  ! interface Quadtree
  !   module procedure init_QuadTree
  ! end interface Quadtree


contains
  !-------------------------------------------------------
  !- constructor for QuadTree
  !-------------------------------------------------------
  type(QuadTree) function init_QuadTree(max_level)
    integer,intent(in) :: max_level
    integer :: num_nodes
    type(Node),allocatable :: nodes(:)
    num_nodes = (4**(max_level+1)-1)/3 
    !print *,"-   initialize QuadTree, num_nodes:",num_nodes
    ! ↓レベルあげるとクッソ重くなる
    allocate(nodes(num_nodes))
    nodes(1)%level = 0
    nodes(:)%is_leaf = .false.
    nodes(:)%num_vor = 0
    init_QuadTree%nodes = nodes
    init_QuadTree%level_hist(:) = 0
  end function init_QuadTree
  !-------------------------------------------------------

  ! Othre subroutine
  function  quadrant_pos(pos,c_pos) result(quad_index)
    double precision,intent(in) :: pos(2) 
    double precision,intent(in) :: c_pos(2)
    integer :: i, quad_index
    quad_index = 0
    if (pos(1) >= c_pos(1)) then
      quad_index = ibset(quad_index,0)
    else 
      quad_index = ibclr(quad_index,0)
    end if
    if (pos(2) >= c_pos(2)) then
      quad_index = ibclr(quad_index,1)
    else 
      quad_index = ibset(quad_index,1)
    end if
    quad_index = quad_index + 1
    !print *,quad_index,pos
  end function quadrant_pos

  subroutine make_tree(vor,qtree,num_p,nterms,theta,t_level,r_min)
    type(Vortex),intent(in) :: vor(:)
    type(QuadTree), intent(out) :: qtree
    integer,intent(in) :: num_p
    integer,intent(in) :: nterms
    double precision,intent(in) :: theta
    integer,intent(in) :: t_level
    double precision,intent(in) :: r_min
    !type(QuadTree),target :: qtree


    ! integer,intent(in) :: max_level
    integer :: i
    integer :: num_nodes
    type(Node),allocatable :: nodes(:)
    integer,save :: max_level = 1

    !print *,"- make_tree"
    ! ----------------------------------
    ! - constructor
    ! ----------------------------------
    qtree%stat = 1 
    do while (qtree%stat == 1)
      qtree%max_level = max_level

      num_nodes = (4**(max_level+1)-1)/3 
      !print *,"-   initialize QuadTree, num_nodes:",num_nodes
      ! ↓レベルあげるとクッソ重くなる
      ! allocate(nodes(num_nodes))
      allocate(qtree%nodes(num_nodes))
      qtree%nodes(1)%level = 0
      qtree%nodes(:)%is_leaf = .false.
      qtree%nodes(:)%num_vor = 0
      ! qtree%nodes = nodes
      qtree%level_hist(:) = 0
      ! ----------------------------------

      ! qtree = QuadTree(max_level)
      qtree%nterms = nterms
      qtree%theta = theta
      qtree%t_level = t_level
      qtree%r_min = r_min 
      !p_qtree => qtree

      qtree%nodes(1)%num_vor = size(vor)
      qtree%nodes(1)%level = 0
      do i = 1,size(qtree%nodes)
        allocate(qtree%nodes(i)%inner_coeffs(qtree%nterms))
        allocate(qtree%nodes(i)%outer_coeffs(qtree%nterms))
        qtree%nodes(i)%outer_coeffs = cmplx(0.d0,0.d0)
        qtree%nodes(i)%inner_coeffs = cmplx(0.d0,0.d0)
      end do

      call add_vortex(vor,qtree)
      !call split_cell(vor,qtree,num_p)
      qtree%stat = 0
      call split_cell(vor,qtree, num_p, 1, .true.)
      if (qtree%stat == -1) then
        deallocate(qtree%nodes)
        qtree%stat = 1
        max_level = max_level + 1
      end if
    end do
    !call interaction_set(p_qtree)
  
  end subroutine make_tree
  
  
  subroutine add_vortex(vor,qtree)
    integer :: i 
    type(QuadTree), intent(inout) :: qtree
    type(Vortex), intent(in) :: vor(:)
    double precision :: region(2)
    double precision :: minValue
    !print *,"-   add_vortex"
    qtree%nodes(1)%num_vor = size(vor)
    allocate(qtree%nodes(1)%vor_id(size(vor)))
    do i = 1,size(vor)
      qtree%nodes(1)%vor_id(i) = i 
    end do
    !+0.1はセル格子とかぶった時にpos2indexが正しく計算されないので、少し広げる
    region(1) = abs(maxval(vor%pos(1)) - minval(vor%pos(1)))+0.0001
    region(2) = abs(maxval(vor%pos(2)) - minval(vor%pos(2)))+0.0001
    if (region(1) > region(2)) then
      qtree%nodes(1)%region = region(1)
    else 
      qtree%nodes(1)%region = region(2) 
    end if
    qtree%nodes(1)%center_pos(1) = 0.5*(maxval(vor(:)%pos(1)) + minval(vor(:)%pos(1)))
    qtree%nodes(1)%center_pos(2) = 0.5*(maxval(vor(:)%pos(2)) + minval(vor(:)%pos(2)))
  end subroutine add_vortex


  subroutine calc_CenterPos(quad_i,p_centerPos, region, c_centerPos)
    implicit none
    integer,intent(in) :: quad_i
    double precision,intent(in) :: p_centerPos(2)
    double precision,intent(in) :: region
    double precision,intent(out) :: c_centerPos(2)
    if (btest(quad_i-1,0)) then
      c_centerPos(1) = p_centerPos(1) + 0.25*region
    else 
      c_centerPos(1) = p_centerPos(1) - 0.25*region
    end if
    if (btest(quad_i-1,1)) then
      c_centerPos(2) = p_centerPos(2) - 0.25*region
    else 
      c_centerPos(2) = p_centerPos(2) + 0.25*region
    end if
  end subroutine calc_centerPos

  
  recursive subroutine split_cell(vor, qtree, num_p, p_index,flag)
    ! cf_ / child_first, p_ / parent
    use m_morton_order
    implicit none
    type(Vortex), intent(in) :: vor(:)
    type(QuadTree), intent(inout),target :: qtree
    integer,intent(in) :: num_p
    integer,intent(in) :: p_index
    logical,intent(in) :: flag
    integer :: c_index(4)
    integer :: i,j,k
    integer :: quad_index 
    integer :: p_num_vor
    integer :: node_index
    integer :: vor_index
    integer,allocatable :: p_quad(:,:)
    integer :: num_p_quad(4)
    type(Vortex) :: tmp_vor
    ! integer,allocatable :: vor_id(:)
    double precision :: x,y
    complex(kind(0d0)) :: z_j,z_star
    double precision,parameter :: pi = 2*acos(0.d0)
    type(Node),pointer :: c_cell


    ! quadrant loop
    if (qtree%nodes(p_index)%level + 1 > qtree%max_level) then
      qtree%stat = -1
      return 
      !print *,"Error: You cannot refer to levels higher than expected."
      !stop
    end if

    p_num_vor = qtree%nodes(p_index)%num_vor
    call getChild(p_index, c_index)
    do i = 1,4 
      allocate(qtree%nodes(c_index(i))%vor_id(p_num_vor))
      !qtree%nodes(c_index(i))%num_vor = 0
    end do

    ! vortex loop
    do i = 1,qtree%nodes(p_index)%num_vor
      tmp_vor = vor(qtree%nodes(p_index)%vor_id(i))
      quad_index = quadrant_pos(tmp_vor%pos,qtree%nodes(p_index)%center_pos)
      node_index = c_index(quad_index)
      qtree%nodes(node_index)%num_vor = qtree%nodes(node_index)%num_vor + 1
      qtree%nodes(node_index)%vor_id(qtree%nodes(node_index)%num_vor) = qtree%nodes(p_index)%vor_id(i)
    end do

    ! quadrant loop
    do i = 1, 4
      call  calc_centerPos(i, &
                           & qtree%nodes(p_index)%center_pos, &
                           & qtree%nodes(p_index)%region, &
                           & qtree%nodes(c_index(i))%center_pos)

      qtree%nodes(c_index(i))%id = c_index(i)
      qtree%nodes(c_index(i))%level = qtree%nodes(p_index)%level + 1
      qtree%nodes(c_index(i))%region = 0.5*qtree%nodes(p_index)%region
      if (qtree%nodes(c_index(i))%num_vor < num_p) then
        if (flag) then
          qtree%nodes(c_index(i))%is_leaf = .true.
          qtree%level_hist(qtree%nodes(c_index(i))%level) = &
              & qtree%level_hist(qtree%nodes(c_index(i))%level) + qtree%nodes(c_index(i))%num_vor
        end if
        if (qtree%nodes(c_index(i))%level < qtree%max_level) then
          ! already reached reaf-cell
          if (qtree%stat /= -1) call split_cell(vor, qtree, num_p, c_index(i), .false.)
        end if
      else
        ! not reach reaf-cell
        if (qtree%stat /= -1) call split_cell(vor, qtree, num_p, c_index(i), .true.)
      end if

    end do
    !print *,"-----------------------"

  end subroutine split_cell

end module m_QuadTree


