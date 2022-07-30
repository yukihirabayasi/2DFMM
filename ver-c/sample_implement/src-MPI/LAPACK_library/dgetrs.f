      subroutine dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

      character          trans
      integer            info, lda, ldb, n, nrhs

      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )

      double precision   one
      parameter          ( one = 1.0d+0 )

      logical            notran

      logical            lsame
      external           lsame

      external           dlaswp, dtrsm, xerbla

      intrinsic          max

      info = 0
      notran = lsame( trans, 'n' )
      if( .not.notran .and. .not.lsame( trans, 't' ) .and. .not.  &
     &    lsame( trans, 'c' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( nrhs.lt.0 ) then
         info = -3
      else if( lda.lt.max( 1, n ) ) then
         info = -5
      else if( ldb.lt.max( 1, n ) ) then
         info = -8
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetrs', -info )
         return
      end if

      if( n.eq.0 .or. nrhs.eq.0 )  &
     &   return

      if( notran ) then

         call dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )

         call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs,  &
     &               one, a, lda, b, ldb )

         call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n,  &
     &               nrhs, one, a, lda, b, ldb )
      else

         call dtrsm( 'left', 'upper', 'transpose', 'non-unit', n, nrhs,  &
     &               one, a, lda, b, ldb )

         call dtrsm( 'left', 'lower', 'transpose', 'unit', n, nrhs, one,  &
     &               a, lda, b, ldb )

         call dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      end if

      return

      end
