      integer function idamax(n,dx,incx)

      integer incx,n

      double precision dx(*)

      double precision dmax
      integer i,ix

      intrinsic dabs

      idamax = 0
      if (n.lt.1 .or. incx.le.0) return
      idamax = 1
      if (n.eq.1) return
      if (incx.eq.1) then

         dmax = dabs(dx(1))
         do i = 2,n
            if (dabs(dx(i)).gt.dmax) then
               idamax = i
               dmax = dabs(dx(i))
            end if
         end do
      else

         ix = 1
         dmax = dabs(dx(1))
         ix = ix + incx
         do i = 2,n
            if (dabs(dx(ix)).gt.dmax) then
               idamax = i
               dmax = dabs(dx(ix))
            end if
            ix = ix + incx
         end do
      end if
      return
      end
