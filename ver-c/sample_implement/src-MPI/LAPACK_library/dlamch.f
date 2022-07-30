      double precision function dlamch( cmach )

      character          cmach

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   rnd, eps, sfmin, small, rmach

      logical            lsame
      external           lsame

      intrinsic          digits, epsilon, huge, maxexponent,  &
     &                   minexponent, radix, tiny

      rnd = one

      if( one.eq.rnd ) then
         eps = epsilon(zero) * 0.5
      else
         eps = epsilon(zero)
      end if

      if( lsame( cmach, 'e' ) ) then
         rmach = eps
      else if( lsame( cmach, 's' ) ) then
         sfmin = tiny(zero)
         small = one / huge(zero)
         if( small.ge.sfmin ) then

            sfmin = small*( one+eps )
         end if
         rmach = sfmin
      else if( lsame( cmach, 'b' ) ) then
         rmach = radix(zero)
      else if( lsame( cmach, 'p' ) ) then
         rmach = eps * radix(zero)
      else if( lsame( cmach, 'n' ) ) then
         rmach = digits(zero)
      else if( lsame( cmach, 'r' ) ) then
         rmach = rnd
      else if( lsame( cmach, 'm' ) ) then
         rmach = minexponent(zero)
      else if( lsame( cmach, 'u' ) ) then
         rmach = tiny(zero)
      else if( lsame( cmach, 'l' ) ) then
         rmach = maxexponent(zero)
      else if( lsame( cmach, 'o' ) ) then
         rmach = huge(zero)
      else
         rmach = zero
      end if

      dlamch = rmach
      return

      end

