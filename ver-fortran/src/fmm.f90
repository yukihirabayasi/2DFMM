module m_fmm
  use m_morton_order
  use m_read_s_file
  !use m_exchange_vortex
  implicit none
contains
  subroutine prepare_fmm(vor,qtree)
    use m_QuadTree
    use m_standard_collection
    implicit none
    type(Vortex),intent(out),allocatable :: vor(:)
    type(QuadTree),intent(out) :: qtree
    integer :: nterms
    integer :: num_p
    double precision :: theta ! for treecode
    integer ::t_level
    double precision :: r_min
    integer :: i

    integer t1, t2, t_rate, t_max, diff
  
    !----------------------------------------------
    !- Read Parameter
    !---------------------------------------------
    open(230,file="../cal_cond/fmm_param.dat")
    read(230,*) nterms
    read(230,*) num_p
    read(230,*) theta
    read(230,*) t_level
    read(230,*) r_min
    close(230)

    call read_s_file(vor)
    !call exchange_vortex(vor)
    call make_tree(vor,qtree,num_p,nterms,theta,t_level,r_min)
    call S2M(qtree,vor)
    call M2M(qtree)
    !call M2L(qtree)
    !call L2L(qtree)

  end subroutine prepare_fmm

  subroutine S2M(qtree,vor)
    use m_QuadTree
    use m_standard_collection
    implicit none
    type(QuadTree),intent(inout),target :: qtree
    type(Vortex), intent(in) :: vor(:)
    complex(kind(0d0)) :: z_j,z_star
    type(Node),pointer :: c_cell
    integer i,j,k

    ! vor_id = qtree%nodes(c_index(i))%vor_id
    do i = 1,size(qtree%nodes)
      if(qtree%nodes(i)%is_leaf) then
        z_star = cmplx(qtree%nodes(i)%center_pos(1),qtree%nodes(i)%center_pos(2))
        c_cell => qtree%nodes(i)
        do k = 0, qtree%nterms-1
          do j = 1, qtree%nodes(i)%num_vor
            ! c_cell%outer_coeffs(k+1) = c_cell%outer_coeffs(k+1) + vor(vor_id(j))%cir*(z_j - z_star)**k
            z_j = cmplx(vor(qtree%nodes(i)%vor_id(j))%pos(1), &
                      & vor(qtree%nodes(i)%vor_id(j))%pos(2))
            c_cell%outer_coeffs(k+1) = c_cell%outer_coeffs(k+1) + &
                      & vor(qtree%nodes(i)%vor_id(j))%cir*(z_j - z_star)**k
          end do
        end do
      end if
    end do
    !------------------------------------------------------------
  
  end subroutine S2M

  
  subroutine M2M(qtree)
    use m_QuadTree
    use m_standard_collection
    implicit none
    type(QuadTree),intent(inout) :: qtree
    integer :: i,j,k,l,n,p_index
    complex(kind(0d0)) :: coeffs(qtree%nterms),p_pos,c_pos
    integer :: level
    double precision :: u_fmm,v_fmm

    do i = size(qtree%nodes),5,-4
      call getParent(i,p_index)
      p_pos = cmplx(qtree%nodes(p_index)%center_pos(1),qtree%nodes(p_index)%center_pos(2))
      do j = 0,-3,-1
        if (qtree%nodes(i+j)%outer_coeffs(1) == cmplx(0.d0,0.d0)) cycle
        c_pos = cmplx(qtree%nodes(i+j)%center_pos(1),qtree%nodes(i+j)%center_pos(2))
        coeffs = cmplx(0.d0,0.d0)
        do n = 0, qtree%nterms-1
          do k = 0, n
            coeffs(n+1) = coeffs(n+1) + qtree%nodes(i+j)%outer_coeffs(k+1)*combinate(n,k)*(-(p_pos-c_pos))**(n-k)
          end do
        end do
        qtree%nodes(p_index)%outer_coeffs = qtree%nodes(p_index)%outer_coeffs + coeffs
      end do
      !if (.not. qtree%nodes(p_index)%is_leaf) qtree%nodes(p_index)%outer_coeffs = coeffs
    end do
    
  end subroutine M2M

  
  subroutine M2L(qtree)
    use m_QuadTree
    use m_standard_collection
    implicit none
    type(QuadTree),intent(inout) :: qtree
    integer :: interaction_set(27)
    integer :: i,j,k,l
    integer :: p_index
    complex(kind(0d0)) :: weight(qtree%nterms,qtree%nterms)
    complex(kind(0d0)) :: c_pos, n_pos, p_pos, coeffs(qtree%nterms),p_coeffs(qtree%nterms)
    complex(kind(0d0)) :: t

    do i = 6, size(qtree%nodes)
      ! from first cell in level 2
      c_pos = cmplx(qtree%nodes(i)%center_pos(1),qtree%nodes(i)%center_pos(2))
      coeffs = cmplx(0.d0,0.d0)
      interaction_set = 0
      call get_interaction_set(i,interaction_set)
      do j = 1, 27
        if(interaction_set(j) == 0) exit
        if(qtree%nodes(interaction_set(j))%outer_coeffs(1) == 0.0d0) cycle
        n_pos = cmplx(qtree%nodes(interaction_set(j))%center_pos(1),qtree%nodes(interaction_set(j))%center_pos(2))

        weight = cmplx(0.d0,0.d0)
        t = c_pos-n_pos
        do k = 0, qtree%nterms-1
          do l = 0, qtree%nterms-1
            weight(l+1,k+1) = ((-1)**k) * combinate(k+l,k)*(t**(-(k+l+1)))
            !weight(l+1,k+1) = ((-1)**k) * combinate(k+l,k)/(t**(-(k+l+1)))
          end do
        end do
        coeffs = coeffs + matmul(qtree%nodes(interaction_set(j))%outer_coeffs,weight) 
      end do
      qtree%nodes(i)%inner_coeffs = coeffs
    end do
  end subroutine M2L


  subroutine L2L(qtree)
    use m_QuadTree
    use m_standard_collection
    implicit none
    type(QuadTree),intent(inout) :: qtree
    integer :: i,j,k,l
    integer :: p_index
    complex(kind(0d0)) :: weight(qtree%nterms,qtree%nterms)
    complex(kind(0d0)) :: c_pos, n_pos, p_pos, coeffs(qtree%nterms),p_coeffs(qtree%nterms)
    complex(kind(0d0)) :: t

    do i = 21, size(qtree%nodes)
      call getParent(i,p_index)
      c_pos = cmplx(qtree%nodes(i)%center_pos(1),qtree%nodes(i)%center_pos(2))
      p_pos = cmplx(qtree%nodes(p_index)%center_pos(1),qtree%nodes(p_index)%center_pos(2))
      weight = cmplx(0.d0,0.d0)
      t = c_pos-p_pos
      do k = 0, qtree%nterms-1
        do l = k, qtree%nterms-1
          weight(l+1,k+1) = combinate(l,k)*(t**(l-k))
        end do
      end do
      p_coeffs = matmul(qtree%nodes(p_index)%inner_coeffs,weight)
      qtree%nodes(i)%inner_coeffs = qtree%nodes(i)%inner_coeffs + p_coeffs
    end do

  end subroutine L2L
  
  subroutine biot_FMM(x,y,qtree,vor,u,v)
    use m_QuadTree
    implicit none
    double precision, intent(in) :: x,y
    type(QuadTree), intent(in) :: qtree
    type(Vortex), intent(in) :: vor(:)
    double precision, intent(out) :: u,v
    integer :: i,j,k
    complex(kind(0d0)) :: phi,z_i,z_star
    double precision :: ganma,first_term
    integer :: t_level,t_index
    integer :: interaction_set(27)
    integer :: neighbor_set(8)

    double precision,parameter :: pi = 2*acos(0.d0)

    double precision :: rx,ry,r,rv2,xai,rv,ax,ay
    double precision :: cell_size
    integer :: t_pos(2)
    double precision :: u_direct(2)

    t_level = qtree%t_level 
    cell_size = qtree%nodes(1)%region / 2**t_level
    t_pos(1) = int((x-qtree%nodes(1)%center_pos(1)+0.5*qtree%nodes(1)%region)/cell_size)
    t_pos(2) = int(-(y-(qtree%nodes(1)%center_pos(2)+0.5*qtree%nodes(1)%region))/cell_size)
    call pos2index(t_pos,t_level,t_index)


    z_i = cmplx(x,y)
    z_star = cmplx(qtree%nodes(t_index)%center_pos(1),qtree%nodes(t_index)%center_pos(2))
    phi = cmplx(0.d0,0.d0)
    do k = 0, qtree%nterms-1
      phi = phi + ((z_i-z_star)**k)*qtree%nodes(t_index)%inner_coeffs(k+1)
    end do
    phi = 0.5d0*phi/pi
    !phi = conjg(phi)
    phi = cmplx(-aimag(phi),real(phi))
    
    u_direct = 0.d0
    call get_neighbor_set(t_index,neighbor_set)
    do i = 1, 8
      if (neighbor_set(i) == 0) cycle
      call direct_evaluate(x,y,vor,neighbor_set(i),qtree,u_direct)
    end do
    !call get_interaction_set(t_index,interaction_set)
    !do i = 1, 27
    !  if (interaction_set(i) == 0) exit
    !  call direct_evaluate(x,y,vor,interaction_set(i),qtree,u_direct)
    !end do
    call direct_evaluate(x,y,vor,t_index,qtree,u_direct)
    u = real(phi) + u_direct(1)
    v = -aimag(phi) + u_direct(2)

 end subroutine biot_FMM
    
        
  subroutine biot_tree(x,y,qtree,vor,u,v) ! evaluate potential
    use m_QuadTree
    double precision,intent(in) :: x,y
    type(QuadTree),intent(in) :: qtree
    type(Vortex),intent(in) :: vor(:)
    double precision,intent(out) :: u,v
    double precision :: u_direct(2)
    complex(kind(0d0)) :: phi
    double precision,parameter :: pi = 2.d0*acos(0.d0)

    phi = cmplx(0.d0,0.d0)
    u_direct = 0.d0
    call tree_evaluate(x, y, vor, 1, qtree, phi, u_direct)
    phi = 0.5d0*cmplx(-aimag(phi),real(phi))/pi
    u = real(phi) + u_direct(1)
    v = -aimag(phi) + u_direct(2)
  end subroutine biot_tree


  recursive subroutine tree_evaluate(x, y, vor, p_index, qtree, phi, u_direct)
    use m_QuadTree
    double precision,intent(in) :: x,y
    type(Vortex),intent(in) :: vor(:)
    integer,intent(in) :: p_index
    type(QuadTree),intent(in),target :: qtree
    complex(kind(0d0)),intent(inout) :: phi
    double precision,intent(inout) :: u_direct(2)
    integer :: c_index(4)
    double precision :: r ! distance between the target point and the center of the cell

    type(Node),pointer :: c_cell
    integer :: i,k
    complex(kind(0d0)) :: z_i,z_star

    z_i = cmplx(x,y)

    if (qtree%nodes(p_index)%is_leaf) then
      ! 指定ノード内の渦要素による影響を計算
      ! reafセルであればvor_idを持っていることを利用
      call direct_evaluate(x,y,vor,p_index,qtree,u_direct)
    else 
      !call testChild(p_index, c_index,qtree%nodes(p_index)%level+1) ! ver don't seek level 1.25s -> 1.24s
      call getChild(p_index, c_index)
      do i = 1, 4
        c_cell => qtree%nodes(c_index(i))
        r = sqrt((x-c_cell%center_pos(1))**2 + (y-c_cell%center_pos(2))**2)
        if (qtree%r_min > r) then
          call direct_evaluate(x,y,vor,c_index(i),qtree,u_direct)
        else
          if (c_cell%region > qtree%theta*r) then
            ! near-field child cell
            call tree_evaluate(x, y, vor, c_index(i), qtree, phi, u_direct)
          else 
            ! far-field child cell (multipole to potential)
            z_star = cmplx(c_cell%center_pos(1),c_cell%center_pos(2))
            do k = 0, qtree%nterms-1
              phi = phi + (z_i - z_star)**(-k-1) * c_cell%outer_coeffs(k+1)
            end do
          end if
        end if
      end do
    end if 

  end subroutine tree_evaluate

  subroutine direct_evaluate(x,y,vor,t_index,qtree,u_direct)
    ! Directly calculate the effects from all the vortices contained in the cell
    use m_QuadTree
    double precision,intent(in) :: x,y
    type(Vortex),intent(in) :: vor(:)
    integer,intent(in) :: t_index ! target index
    type(QuadTree),intent(in) :: qtree
    double precision,intent(inout) :: u_direct(2)
    integer :: i
    
    double precision,parameter :: pi = 2.d0*acos(0.d0)
    double precision :: rx,ry,r,rv2,xai,rv,ax,ay
    ax = 0.d0
    ay = 0.d0
    do i = 1, qtree%nodes(t_index)%num_vor
      rx = x - vor(qtree%nodes(t_index)%vor_id(i))%pos(1)
      ry = y - vor(qtree%nodes(t_index)%vor_id(i))%pos(2)
      r = sqrt(rx**2 + ry**2)
      if ( r > 1.0d-10) then
        rv2 = 1.0d0 / r**2
        xai = r / vor(qtree%nodes(t_index)%vor_id(i))%cor
        rv  = vor(qtree%nodes(t_index)%vor_id(i))%cir * ( 1.0d0 - exp( -xai*xai ) )
        ax = ax + rv*ry*rv2
        ay = ay - rv*rx*rv2
      end if
    end do
    u_direct(1) = u_direct(1) + 0.5d0*ax/pi
    u_direct(2) = u_direct(2) + 0.5d0*ay/pi
  end subroutine direct_evaluate
end module m_fmm
