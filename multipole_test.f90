module m_Vortex
  implicit none
  type :: Vortex
    double precision :: pos(2)
    double precision :: cir
  end type Vortex
end module m_Vortex

program test
  use m_Vortex
  implicit none
  double precision :: u1,v1
  double precision :: u2,v2
  double precision :: x,y
  type(Vortex) :: vor

  vor%pos(1) = 1.0d0
  vor%pos(2) = 1.0d0
  vor%cir = 1.0d0

  x = 3.0d0
  y = 3.0d0

  call biot_b(x,y,u1,v1,vor)
  call multipole(x,y,vor,u2,v2)

  print *,"biot_b:",u1,v1
  print *,"multipole:",u2,v2

contains
  subroutine biot_b(x,y,u1,v1,vor)
    implicit none
    double precision,intent(in) :: x,y
    double precision,intent(out) :: u1,v1
    type(Vortex),intent(in) :: vor

    double precision :: r,rx,ry,rv2,rv
    
    rx = x - vor%pos(1)
    ry = y - vor%pos(2)
    r = sqrt(rx**2 + ry**2)
    if ( r > 1.0d-6) then
      rv2 = 1.0d0 / r**2
      !xai = r / vor%co
      rv = vor%cir
      u1 = u1 + ry*rv*rv2
      v1 = v1 - rx*rv*rv2
    end if
  end subroutine biot_b

  subroutine multipole(x,y,vor,u2,v2)
    implicit none
    double precision,intent(in) :: x,y
    type(Vortex) :: vor
    double precision,intent(out) :: u2,v2
    
    complex(kind(0d0)) :: phi
    complex(kind(0d0)) :: z,zi
    complex(kind(0d0)) :: temp
    complex(kind(0d0)) :: sec_term = (0d0,0d0)
    double precision :: q
    double precision,parameter :: pi = 2.0d0*acos(0.0d0)
    integer ::i, k = 6

    z = cmplx(x,y)
    zi = cmplx(vor%pos(1),vor%pos(2))
    q = cmplx(0.d0,-0.5*vor%cir/pi)
    do i = 1,k
      temp = (0.5d0*vor%cir*zi**i)/(i*pi*z**i)
      print *,z
      temp = cmplx(-aimag(temp),real(temp))
      sec_term = sec_term + temp
    end do

    phi = q*log(z) + sec_term

    u2 = real(phi)
    v2 = -aimag(phi)

  end subroutine multipole
end program test
