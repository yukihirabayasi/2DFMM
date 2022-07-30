      subroutine dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )

      integer            info, lda, ldb, n, nrhs

      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )

      external           dgetrf, dgetrs, xerbla

      intrinsic          max

      info = 0
      if( n.lt.0 ) then
         info = -1
      else if( nrhs.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, n ) ) then
         info = -4
      else if( ldb.lt.max( 1, n ) ) then
         info = -7
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgesv ', -info )
         return
      end if

      call dgetrf( n, n, a, lda, ipiv, info )
      if( info.eq.0 ) then

         call dgetrs( 'no transpose', n, nrhs, a, lda, ipiv, b, ldb,  &
     &                info )
      end if
      return

      end
