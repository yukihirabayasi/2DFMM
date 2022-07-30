      subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)

      double precision alpha
      integer incx,incy,lda,m,n

      double precision a(lda,*),x(*),y(*)

      double precision zero
      parameter (zero=0.0d+0)

      double precision temp
      integer i,info,ix,j,jy,kx

      external xerbla

      intrinsic max

      info = 0
      if (m.lt.0) then
          info = 1
      else if (n.lt.0) then
          info = 2
      else if (incx.eq.0) then
          info = 5
      else if (incy.eq.0) then
          info = 7
      else if (lda.lt.max(1,m)) then
          info = 9
      end if
      if (info.ne.0) then
          call xerbla('dger  ',info)
          return
      end if

      if ((m.eq.0) .or. (n.eq.0) .or. (alpha.eq.zero)) return

      if (incy.gt.0) then
          jy = 1
      else
          jy = 1 - (n-1)*incy
      end if
      if (incx.eq.1) then
          do 20 j = 1,n
              if (y(jy).ne.zero) then
                  temp = alpha*y(jy)
                  do 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
   10             continue
              end if
              jy = jy + incy
   20     continue
      else
          if (incx.gt.0) then
              kx = 1
          else
              kx = 1 - (m-1)*incx
          end if
          do 40 j = 1,n
              if (y(jy).ne.zero) then
                  temp = alpha*y(jy)
                  ix = kx
                  do 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
   30             continue
              end if
              jy = jy + incy
   40     continue
      end if

      return

      end
