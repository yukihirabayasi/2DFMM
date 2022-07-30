module m_tree_bindc
  use,intrinsic :: iso_c_binding 
  interface 
    subroutine make_tree(nvor_b,vor_b,cir_b,cor_b) bind(c)
      import
      integer(c_int),intent(in),value :: nvor_b
      real(c_double),intent(in) :: vor_b(2,nvor_b),cir_b(nvor_b),cor_b(nvor_b)
    end subroutine

    subroutine biot_tree(x, y, u, v) bind(c)
      import
      real(c_double),intent(in),value  :: x,y
      real(c_double),intent(out) :: u,v
    end subroutine

  end interface
end module

