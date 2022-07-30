include './includes/modules.f90'
include './read_s_file.f90'
include './m_tree_bindc.f90'
include './includes/d_biot.f90'

program main
  use global
  use m_tree_bindc
  implicit none

  real(8) :: x,y
  real(8) :: u,v
  real(8) :: ax,ay

  integer nterms
  real(8) theta_tree
  real(8) r_min

  x = 0.72d0
  y = -0.01d0

  nterms = 12
  theta_tree = 0.4
  r_min = 0.01

  call read_s_file

  
  call make_tree(nvor_b, &
               & vor_b(:,1:nvor_b), &
               & cir_b(1:nvor_b),   &
               & cor_b(1:nvor_b))   
  call biot_tree(x,y,u,v)
  print *,"u,v",u,v

  ax = 0
  ay = 0
  call biot_b(x,y,ax,ay)
  print *,"ax,ay",ax,ay
end program

