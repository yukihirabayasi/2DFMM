      integer          function ieeeck( ispec, zero, one )

      integer            ispec
      real               one, zero

      real               nan1, nan2, nan3, nan4, nan5, nan6, neginf,  &
     &                   negzro, newzro, posinf

      ieeeck = 1

      posinf = one / zero
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if

      neginf = -one / zero
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if

      negzro = one / ( neginf+one )
      if( negzro.ne.zero ) then
         ieeeck = 0
         return
      end if

      neginf = one / negzro
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if

      newzro = negzro + zero
      if( newzro.ne.zero ) then
         ieeeck = 0
         return
      end if

      posinf = one / newzro
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if

      neginf = neginf*posinf
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if

      posinf = posinf*posinf
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if

      if( ispec.eq.0 )  &
     &   return

      nan1 = posinf + neginf

      nan2 = posinf / neginf

      nan3 = posinf / posinf

      nan4 = posinf*zero

      nan5 = neginf*negzro

      nan6 = nan5*zero

      if( nan1.eq.nan1 ) then
         ieeeck = 0
         return
      end if

      if( nan2.eq.nan2 ) then
         ieeeck = 0
         return
      end if

      if( nan3.eq.nan3 ) then
         ieeeck = 0
         return
      end if

      if( nan4.eq.nan4 ) then
         ieeeck = 0
         return
      end if

      if( nan5.eq.nan5 ) then
         ieeeck = 0
         return
      end if

      if( nan6.eq.nan6 ) then
         ieeeck = 0
         return
      end if

      return
      end
