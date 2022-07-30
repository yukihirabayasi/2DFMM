      subroutine xerbla( srname, info )

      character*(*)      srname
      integer            info

      intrinsic          len_trim

      write( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info

      stop

 9999 format( ' ** on entry to ', a, ' parameter number ', i2, ' had ',  &
     &      'an illegal value' )

      end
