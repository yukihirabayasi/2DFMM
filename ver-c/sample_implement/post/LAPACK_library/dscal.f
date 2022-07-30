      subroutine dscal(n,da,dx,incx)

      double precision da
      integer incx,n

      double precision dx(*)

      integer i,m,mp1,nincx

      intrinsic mod

      if (n.le.0 .or. incx.le.0) return
      if (incx.eq.1) then

         m = mod(n,5)
         if (m.ne.0) then
            do i = 1,m
               dx(i) = da*dx(i)
            end do
            if (n.lt.5) return
         end if
         mp1 = m + 1
         do i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         end do
      else

         nincx = n*incx
         do i = 1,nincx,incx
            dx(i) = da*dx(i)
         end do
      end if
      return
      end
