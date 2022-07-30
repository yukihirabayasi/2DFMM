      subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )

      integer            incx, k1, k2, lda, n

      integer            ipiv( * )
      double precision   a( lda, * )

      integer            i, i1, i2, inc, ip, ix, ix0, j, k, n32
      double precision   temp

      if( incx.gt.0 ) then
         ix0 = k1
         i1 = k1
         i2 = k2
         inc = 1
      else if( incx.lt.0 ) then
         ix0 = k1 + ( k1-k2 )*incx
         i1 = k2
         i2 = k1
         inc = -1
      else
         return
      end if

      n32 = ( n / 32 )*32
      if( n32.ne.0 ) then
         do 30 j = 1, n32, 32
            ix = ix0
            do 20 i = i1, i2, inc
               ip = ipiv( ix )
               if( ip.ne.i ) then
                  do 10 k = j, j + 31
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
   10             continue
               end if
               ix = ix + incx
   20       continue
   30    continue
      end if
      if( n32.ne.n ) then
         n32 = n32 + 1
         ix = ix0
         do 50 i = i1, i2, inc
            ip = ipiv( ix )
            if( ip.ne.i ) then
               do 40 k = n32, n
                  temp = a( i, k )
                  a( i, k ) = a( ip, k )
                  a( ip, k ) = temp
   40          continue
            end if
            ix = ix + incx
   50    continue
      end if

      return

      end
