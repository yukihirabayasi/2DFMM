include './morton_order.f90'
include './quadtree.f90'
include './read_s_file.f90'
include './standard_collection.f90'
include './fmm.f90'
!include './includes/d_biot.f90'


program testFMM
  !use morton_order
  use m_fmm
  use m_QuadTree


  implicit none
  type(QuadTree) :: qtree
  type(Vortex),allocatable :: vor(:)
  double precision :: x,y
  double precision :: u_tree,v_tree
  double precision :: u_fmm,v_fmm
  double precision :: u_direct(2)
  integer t1, t2, t_rate, t_max, diff

  integer :: i

  x = 0.72d0 
  y = -0.01d0


  
  !----------------------------------------------------
  !- Main
  !----------------------------------------------------
  print *,"-- Process Log --------------------"
  call prepare_fmm(vor,qtree)
  print *,"-----------------------------------"

  print *,""
  print *,"-- Result -------------------------"
  call system_clock(t1)
  do i = 1,1!0000
    call biot_FMM(x,y,qtree,vor,u_fmm,v_fmm)
  end do
  call system_clock(t2, t_rate, t_max)  

  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  !call biot_b_M2MTest(qtree)
  print *,"- Calculation(FMM):",u_fmm,v_fmm
  print *,"- Time(FMM):",diff/dble(t_rate)
  print *,""

  call system_clock(t1)
  do i = 1,1!0000
    call biot_tree(x,y,qtree,vor,u_tree,v_tree)
  end do
  call system_clock(t2, t_rate, t_max)  

  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  !call biot_b_M2MTest(qtree)
  !print *,"Calulation using FMM:", u_fmm,v_fmm
  print *,"- Calculation(Tree):",u_tree,v_tree
  print *,"- Time(Tree):",diff/dble(t_rate)
  print *,""

  call system_clock(t1)
  u_direct(:) = 0.d0
  do i = 87382,349525
    call direct_evaluate(x,y,vor,i,qtree,u_direct)
  end do
  call system_clock(t2, t_rate, t_max)  
  if ( t2 < t1 ) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  endif
  print *,"- Calculation(Direct):", u_direct(1),u_direct(2)
  print *,"- Time(Direct):",10*diff/dble(t_rate)
  !print *,""
  !print *,"Error:",abs((u_tree-u_direct)/u_direct)*100.d0,"%",abs((v_tree-v_direct)/v_direct)*100.d0,"%"
  print *,"-----------------------------------"





contains


! ************************************************************************** 
  subroutine biot_b_M2MTest ( qtree )
! ************************************************************************** 
! -------------------------------------------------------------------------- 
! THIS SUBROUTINE IS TO BE USED ONLY UNDER THE STIPULATIONS OF THE LICENSING 
! AGREEMENT.                                                                 
! -------------------------------------------------------------------------- 
!============================================================================
!  Calculate induced velocity at a point(x,y) from a VORTEX BLOB ELEMENT.    
!============================================================================

  use global
  use parameter
  use math_constant
  use m_morton_order
  implicit none
  include './includes/inc_2d'
  type(QuadTree) :: qtree
  complex(kind(0d0)) :: phi,z_i,z_star
  integer :: index = 3297
  double precision :: vx,vy,minX,maxX,minY,maxY

  print *,"--------------------------------"
  print *,"- biot_b_M2M_Test"
  print *,"--------------------------------"

  !print *,qtree%nodes(index)%
  minX = qtree%nodes(index)%center_pos(1)-0.5*qtree%nodes(index)%region
  maxX = qtree%nodes(index)%center_pos(1)+0.5*qtree%nodes(index)%region
  minY = qtree%nodes(index)%center_pos(2)-0.5*qtree%nodes(index)%region
  maxY = qtree%nodes(index)%center_pos(2)+0.5*qtree%nodes(index)%region

  ! Direct Summation
  vx = 0.0d0
  vy = 0.0d0

  xi = 0.d0
  yi = 0.d0

  do i = 1,nvor_b
    if (minX<=vor_b(1,i) .and. vor_b(1,i)<=maxX .and. minY<=vor_b(2,i) .and. vor_b(2,i)<=MaxY) then
      rx = xi - vor_b(1,i)
      ry = yi - vor_b(2,i)
      r  = dsqrt( rx**2 + ry**2 )
      if ( r .gt. 1.0d-6 ) then
        rv2 = 1.0d0 / r**2
        xai = r / cor_b(i)
        rv  = cir_b(i)! * ( 1.0d0 - dexp( -xai**2 ) )
        vx  = vx + ry*rv*rv2
        vy  = vy - rx*rv*rv2
      end if
    end if
  end do
  vx = vx * pitwo
  vy = vy * pitwo
  print *,"- Direct Summation:",vx,vy

  ! FMM Calculation
  z_i = cmplx(xi,yi)
  z_star = cmplx(qtree%nodes(index)%center_pos(1),qtree%nodes(index)%center_pos(2))
  phi = cmplx(0.d0,0.d0)
  do k = 0,qtree%nterms-1
    phi = phi + (z_i - z_star)**(-k-1) * qtree%nodes(index)%outer_coeffs(k+1)
  end do
  phi = 0.5d0*cmplx(-aimag(phi),real(phi))/pi
  print *,"- FMM",real(phi),-aimag(phi)

  print *,"--------------------------------"
  return

  end subroutine biot_b_M2MTest
end program testFMM
