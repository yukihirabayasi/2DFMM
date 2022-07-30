      integer function iparmq( ispec, name, opts, n, ilo, ihi, lwork )

      integer            ihi, ilo, ispec, lwork, n
      character          name*( * ), opts*( * )

      integer            inmin, inwin, inibl, ishfts, iacc22
      parameter          ( inmin = 12, inwin = 13, inibl = 14,  &
     &                   ishfts = 15, iacc22 = 16 )
      integer            nmin, k22min, kacmin, nibble, knwswp
      parameter          ( nmin = 75, k22min = 14, kacmin = 14,  &
     &                   nibble = 14, knwswp = 500 )
      real               two
      parameter          ( two = 2.0 )

      integer            nh, ns
      integer            i, ic, iz
      character          subnam*6

      intrinsic          log, max, mod, nint, real

      if( ( ispec.eq.ishfts ) .or. ( ispec.eq.inwin ) .or.  &
     &    ( ispec.eq.iacc22 ) ) then

         nh = ihi - ilo + 1
         ns = 2
         if( nh.ge.30 )    &
     &      ns = 4
         if( nh.ge.60 )    &
     &      ns = 10
         if( nh.ge.150 )   &
     &      ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
         if( nh.ge.590 )   &
     &      ns = 64
         if( nh.ge.3000 )  &
     &      ns = 128
         if( nh.ge.6000 )  &
     &      ns = 256
         ns = max( 2, ns-mod( ns, 2 ) )
      end if

      if( ispec.eq.inmin ) then

         iparmq = nmin

      else if( ispec.eq.inibl ) then

         iparmq = nibble

      else if( ispec.eq.ishfts ) then

         iparmq = ns

      else if( ispec.eq.inwin ) then

         if( nh.le.knwswp ) then
            iparmq = ns
         else
            iparmq = 3*ns / 2
         end if

      else if( ispec.eq.iacc22 ) then

         iparmq = 0
         subnam = name
         ic = ichar( subnam( 1: 1 ) )
         iz = ichar( 'z' )
         if( iz.eq.90 .or. iz.eq.122 ) then

            if( ic.ge.97 .and. ic.le.122 ) then
               subnam( 1: 1 ) = char( ic-32 )
               do i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  if( ic.ge.97 .and. ic.le.122 )  &
     &               subnam( i: i ) = char( ic-32 )
               end do
            end if

         else if( iz.eq.233 .or. iz.eq.169 ) then

            if( ( ic.ge.129 .and. ic.le.137 ) .or.  &
     &          ( ic.ge.145 .and. ic.le.153 ) .or.  &
     &          ( ic.ge.162 .and. ic.le.169 ) ) then
               subnam( 1: 1 ) = char( ic+64 )
               do i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  if( ( ic.ge.129 .and. ic.le.137 ) .or.         &
     &                ( ic.ge.145 .and. ic.le.153 ) .or.         &
     &                ( ic.ge.162 .and. ic.le.169 ) )subnam( i:  &
     &                i ) = char( ic+64 )
               end do
            end if

         else if( iz.eq.218 .or. iz.eq.250 ) then

            if( ic.ge.225 .and. ic.le.250 ) then
               subnam( 1: 1 ) = char( ic-32 )
               do i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  if( ic.ge.225 .and. ic.le.250 )  &
     &               subnam( i: i ) = char( ic-32 )
               end do
            end if
         end if

         if( subnam( 2:6 ).eq.'gghrd' .or.  &
     &       subnam( 2:6 ).eq.'gghd3' ) then
            iparmq = 1
            if( nh.ge.k22min )  &
     &         iparmq = 2
         else if ( subnam( 4:6 ).eq.'exc' ) then
            if( nh.ge.kacmin )  &
     &         iparmq = 1
            if( nh.ge.k22min )  &
     &         iparmq = 2
         else if ( subnam( 2:6 ).eq.'hseqr' .or.  &
     &             subnam( 2:5 ).eq.'laqr' ) then
            if( ns.ge.kacmin )  &
     &         iparmq = 1
            if( ns.ge.k22min )  &
     &         iparmq = 2
         end if

      else

         iparmq = -1

      end if

      end
