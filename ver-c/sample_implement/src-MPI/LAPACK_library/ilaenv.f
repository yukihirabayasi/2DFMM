      integer function ilaenv( ispec, name, opts, n1, n2, n3, n4 )

      character*( * )    name, opts
      integer            ispec, n1, n2, n3, n4

      integer            i, ic, iz, nb, nbmin, nx
      logical            cname, sname, twostage
      character          c1*1, c2*2, c4*2, c3*3, subnam*16

      intrinsic          char, ichar, int, min, real

      integer            ieeeck, iparmq, iparam2stage
      external           ieeeck, iparmq, iparam2stage

      go to ( 10, 10, 10, 80, 90, 100, 110, 120,  &
     &        130, 140, 150, 160, 160, 160, 160, 160)ispec

      ilaenv = -1
      return

   10 continue

      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1: 1 ) )
      iz = ichar( 'z' )
      if( iz.eq.90 .or. iz.eq.122 ) then

         if( ic.ge.97 .and. ic.le.122 ) then
            subnam( 1: 1 ) = char( ic-32 )
            do 20 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ic.ge.97 .and. ic.le.122 )  &
     &            subnam( i: i ) = char( ic-32 )
   20       continue
         end if

      else if( iz.eq.233 .or. iz.eq.169 ) then

         if( ( ic.ge.129 .and. ic.le.137 ) .or.  &
     &       ( ic.ge.145 .and. ic.le.153 ) .or.  &
     &       ( ic.ge.162 .and. ic.le.169 ) ) then
            subnam( 1: 1 ) = char( ic+64 )
            do 30 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ( ic.ge.129 .and. ic.le.137 ) .or.         &
     &             ( ic.ge.145 .and. ic.le.153 ) .or.         &
     &             ( ic.ge.162 .and. ic.le.169 ) )subnam( i:  &
     &             i ) = char( ic+64 )
   30       continue
         end if

      else if( iz.eq.218 .or. iz.eq.250 ) then

         if( ic.ge.225 .and. ic.le.250 ) then
            subnam( 1: 1 ) = char( ic-32 )
            do 40 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ic.ge.225 .and. ic.le.250 )  &
     &            subnam( i: i ) = char( ic-32 )
   40       continue
         end if
      end if

      c1 = subnam( 1: 1 )
      sname = c1.eq.'s' .or. c1.eq.'d'
      cname = c1.eq.'c' .or. c1.eq.'z'
      if( .not.( cname .or. sname ) )  &
     &   return
      c2 = subnam( 2: 3 )
      c3 = subnam( 4: 6 )
      c4 = c3( 2: 3 )
      twostage = len( subnam ).ge.11  &
     &           .and. subnam( 11: 11 ).eq.'2'

      go to ( 50, 60, 70 )ispec

   50 continue

      nb = 1

      if( c2.eq.'ge' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or.  &
     &            c3.eq.'qlf' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'qr ') then
            if( n3 .eq. 1) then
               if( sname ) then

                  if ((n1*n2.le.131072).or.(n1.le.8192)) then
                     nb = n1
                  else
                     nb = 32768/n2
                  end if
               else
                  if ((n1*n2.le.131072).or.(n1.le.8192)) then
                     nb = n1
                  else
                     nb = 32768/n2
                  end if
               end if
            else
               if( sname ) then
                  nb = 1
               else
                  nb = 1
               end if
            end if
         else if( c3.eq.'lq ') then
            if( n3 .eq. 2) then
               if( sname ) then

                  if ((n1*n2.le.131072).or.(n1.le.8192)) then
                     nb = n1
                  else
                     nb = 32768/n2
                  end if
               else
                  if ((n1*n2.le.131072).or.(n1.le.8192)) then
                     nb = n1
                  else
                     nb = 32768/n2
                  end if
               end if
            else
               if( sname ) then
                  nb = 1
               else
                  nb = 1
               end if
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'po' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( twostage ) then
                  nb = 192
               else
                  nb = 64
               end if
            else
               if( twostage ) then
                  nb = 192
               else
                  nb = 64
               end if
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nb = 32
         else if( sname .and. c3.eq.'gst' ) then
            nb = 64
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trf' ) then
            if( twostage ) then
               nb = 192
            else
               nb = 64
            end if
         else if( c3.eq.'trd' ) then
            nb = 32
         else if( c3.eq.'gst' ) then
            nb = 64
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nb = 32
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nb = 32
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nb = 32
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nb = 32
            end if
         end if
      else if( c2.eq.'gb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'pb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'tr' ) then
         if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if ( c3.eq.'evc' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'la' ) then
         if( c3.eq.'uum' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( sname .and. c2.eq.'st' ) then
         if( c3.eq.'ebz' ) then
            nb = 1
         end if
      else if( c2.eq.'gg' ) then
         nb = 32
         if( c3.eq.'hd3' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         end if
      end if
      ilaenv = nb
      return

   60 continue

      nbmin = 2
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or. c3.eq.  &
     &       'qlf' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nbmin = 8
            else
               nbmin = 8
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nbmin = 2
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nbmin = 2
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nbmin = 2
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nbmin = 2
            end if
         end if
      else if( c2.eq.'gg' ) then
         nbmin = 2
         if( c3.eq.'hd3' ) then
            nbmin = 2
         end if
      end if
      ilaenv = nbmin
      return

   70 continue

      nx = 0
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or. c3.eq.  &
     &       'qlf' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( sname .and. c3.eq.'trd' ) then
            nx = 32
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nx = 32
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nx = 128
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.  &
     &          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )  &
     &           then
               nx = 128
            end if
         end if
      else if( c2.eq.'gg' ) then
         nx = 128
         if( c3.eq.'hd3' ) then
            nx = 128
         end if
      end if
      ilaenv = nx
      return

   80 continue

      ilaenv = 6
      return

   90 continue

      ilaenv = 2
      return

  100 continue

      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
      return

  110 continue

      ilaenv = 1
      return

  120 continue

      ilaenv = 50
      return

  130 continue

      ilaenv = 25
      return

  140 continue

      ilaenv = 1
      if( ilaenv.eq.1 ) then
         ilaenv = ieeeck( 1, 0.0, 1.0 )
      end if
      return

  150 continue

      ilaenv = 1
      if( ilaenv.eq.1 ) then
         ilaenv = ieeeck( 0, 0.0, 1.0 )
      end if
      return

  160 continue

      ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
      return

      end
