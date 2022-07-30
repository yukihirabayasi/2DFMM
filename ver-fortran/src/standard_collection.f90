module m_standard_collection
  interface
     function combinate(n,r) result(ans)
      integer,intent(in) :: n,r
      integer(8) :: ans
    end function combinate 

    function fractorial(a) result(ans)
      integer,intent(in) :: a
      integer(8) :: ans
    end function fractorial
  end interface

end module m_standard_collection

function combinate(n,r) result(ans)
! ans = nCr = n! / (r!(n-r)!)
  integer,intent(in) :: n,r
  integer(8) :: ans
  interface
    function fractorial(a) result(ans)
      integer,intent(in) :: a
      integer(8) :: ans
    end function fractorial
  end interface
  ans = fractorial(n)/(fractorial(r)*fractorial(n-r))
end function combinate

function fractorial(a) result(ans)
! ans = a! = a * (a-1) * ... * 1
  integer,intent(in) :: a
  integer(8) :: i, ans
  if (a > 20) then
    print *,"Error: An overflow has occurred. Factorial argument must be 12 or less."
    stop
  end if
  ans = 1
  do i = 1, a
    ans = ans * i
  end do
end function fractorial
