program main
  use,intrinsic :: iso_c_binding 
  implicit none

  integer,parameter :: n_vor_b = 3000001

  !! fortranは参照渡し、Cは値渡しなので、その点に注意すること。
  !! また、Treeの情報を示す変数は全てグルーバル変数とすることで、fortran側への変数の伝達を最小限に抑えた。
  !! その弊害として、常に同じメモリに格納されたデータを用いるので、
  !! きちんと新しい流れ場を計算する際には、ツリーの情報を初期化する必要がある。

  interface 
    subroutine make_tree(nvor_b,poi_x,poi_y,cir,cor) bind(c)
      import
      integer(c_int),intent(in) :: nvor_b 
      real(c_double),intent(in) :: poi_x(nvor_b),poi_y(nvor_b),cir(nvor_b),cor(nvor_b)
    end subroutine

    subroutine biot_tree(x, y, u, v) bind(c)
      import
      real(c_double),intent(in)  :: x,y
      real(c_double),intent(out) :: u,v
    end subroutine
  end interface

  integer i
  real(8) :: x,y
  real(8) :: u,v
  integer :: nvor_b
  real(8) :: poi_x(n_vor_b), poi_y(n_vor_b)
  real(8) :: cir(n_vor_b), cor(n_vor_b)

  x = 3.0
  y = 2.0

  nvor_b = 10

  do i = 1,nvor_b
    poi_x(i) = real(i)
  end do

  call make_tree(nvor_b, poi_x(1:nvor_b), poi_y(1:nvor_b), cir(1:nvor_b), cor(1:nvor_b))
  call biot_tree(x,y,u,v)
  print *,"u,v",u,v
end program

