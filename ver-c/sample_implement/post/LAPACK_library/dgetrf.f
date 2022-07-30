      subroutine dgetrf ( m, n, a, lda, ipiv, info)

      integer            info, lda, m, n

      integer            ipiv( * )
      double precision   a( lda, * )

      double precision   one
      parameter          ( one = 1.0d+0 )

      integer            i, iinfo, j, jb, nb

      external           dgemm, dgetf2, dlaswp, dtrsm, xerbla

      integer            ilaenv
      external           ilaenv

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
         call xerbla( 'dgetrf', -info )
         return
      end if

      if( m.eq.0 .or. n.eq.0 )  &
     &   return

      nb = ilaenv( 1, 'dgetrf', ' ', m, n, -1, -1 )
      if( nb.le.1 .or. nb.ge.min( m, n ) ) then

         call dgetf2( m, n, a, lda, ipiv, info )
      else

         do 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )

            call dgemm( 'no transpose', 'no transpose',      &
     &                 m-j+1, jb, j-1, -one,                 &
     &                 a( j, 1 ), lda, a( 1, j ), lda, one,  &
     &                 a( j, j ), lda )

            call dgetf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )

            if( info.eq.0 .and. iinfo.gt.0 )  &
     &         info = iinfo + j - 1
            do 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       continue

            call dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )

            if ( j+jb.le.n ) then

               call dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,  &
     &                     ipiv, 1 )

               call dgemm( 'no transpose', 'no transpose',         &
     &                    jb, n-j-jb+1, j-1, -one,                 &
     &                    a( j, 1 ), lda, a( 1, j+jb ), lda, one,  &
     &                    a( j, j+jb ), lda )

               call dtrsm( 'left', 'lower', 'no transpose', 'unit',  &
     &                    jb, n-j-jb+1, one, a( j, j ), lda,         &
     &                    a( j, j+jb ), lda )
            end if

   20    continue

      end if
      return

      end
