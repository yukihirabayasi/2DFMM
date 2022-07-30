!include './includes/modules.f90'

module m_exchange_vortex
  implicit none
contains
  subroutine exchange_vortex(vor)
    use m_QuadTree
  
    use global
    use parameter
    use math_constant
  
    type(Vortex),intent(out),allocatable :: vor(:)
  
    character(50), save :: f_inp
    character(50), save :: f_debug
    integer nf_inp
    integer nf_debug
    integer i,j

    allocate(vor(nvor_b))
    do i = 1, nvor_b
      vor(i)%pos(1) = vor_b(1,i)
      vor(i)%pos(2) = vor_b(2,i)
      vor(i)%cir = cir_b(i)
      vor(i)%cor = cor_b(i)
    end do

  end subroutine exchange_vortex
end module
