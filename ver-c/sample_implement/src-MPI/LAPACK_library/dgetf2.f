      subroutine dgetf2( m, n, a, lda, ipiv, info )

      integer            info, lda, m, n

      integer            ipiv( * )
      double precision   a( lda, * )

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   sfmin
      integer            i, j, jp

      double precision   dlamch
      integer            idamax
      external           dlamch, idamax

      external           dger, dscal, dswap, xerbla

      intrinsic          max, min

      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetf2', -info )
         return
      end if

      if( m.eq.0 .or. n.eq.0 )  &
     &   return

      sfmin = dlamch('s')

      do 10 j = 1, min( m, n )

         jp = j - 1 + idamax( m-j+1, a( j, j ), 1 )
         ipiv( j ) = jp
         if( a( jp, j ).ne.zero ) then

            if( jp.ne.j )  &
     &         call dswap( n, a( j, 1 ), lda, a( jp, 1 ), lda )

            if( j.lt.m ) then
               if( abs(a( j, j )) .ge. sfmin ) then
                  call dscal( m-j, one / a( j, j ), a( j+1, j ), 1 )
               else
                 do 20 i = 1, m-j
                    a( j+i, j ) = a( j+i, j ) / a( j, j )
   20            continue
               end if
            end if

         else if( info.eq.0 ) then

            info = j
         end if

         if( j.lt.min( m, n ) ) then

            call dger( m-j, n-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda,  &
     &                 a( j+1, j+1 ), lda )
         end if
   10 continue
      return

      end
