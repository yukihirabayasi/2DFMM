      subroutine dswap(n,dx,incx,dy,incy)

      integer incx,incy,n

      double precision dx(*),dy(*)

      double precision dtemp
      integer i,ix,iy,m,mp1

      intrinsic mod

      if (n.le.0) return
      if (incx.eq.1 .and. incy.eq.1) then

         m = mod(n,3)
         if (m.ne.0) then
            do i = 1,m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            end do
            if (n.lt.3) return
         end if
         mp1 = m + 1
         do i = mp1,n,3
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i+1)
            dx(i+1) = dy(i+1)
            dy(i+1) = dtemp
            dtemp = dx(i+2)
            dx(i+2) = dy(i+2)
            dy(i+2) = dtemp
         end do
      else

         ix = 1
         iy = 1
         if (incx.lt.0) ix = (-n+1)*incx + 1
         if (incy.lt.0) iy = (-n+1)*incy + 1
         do i = 1,n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         end do
      end if
      return
      end
