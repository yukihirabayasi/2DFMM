      logical function lsame(ca,cb)

      character ca,cb

      intrinsic ichar

      integer inta,intb,zcode

      lsame = ca .eq. cb
      if (lsame) return

      zcode = ichar('z')

      inta = ichar(ca)
      intb = ichar(cb)

      if (zcode.eq.90 .or. zcode.eq.122) then

          if (inta.ge.97 .and. inta.le.122) inta = inta - 32
          if (intb.ge.97 .and. intb.le.122) intb = intb - 32

      else if (zcode.eq.233 .or. zcode.eq.169) then

          if (inta.ge.129 .and. inta.le.137 .or.  &
   &          inta.ge.145 .and. inta.le.153 .or.  &
   &          inta.ge.162 .and. inta.le.169) inta = inta + 64
          if (intb.ge.129 .and. intb.le.137 .or.  &
   &          intb.ge.145 .and. intb.le.153 .or.  &
   &          intb.ge.162 .and. intb.le.169) intb = intb + 64

      else if (zcode.eq.218 .or. zcode.eq.250) then

          if (inta.ge.225 .and. inta.le.250) inta = inta - 32
          if (intb.ge.225 .and. intb.le.250) intb = intb - 32
      end if
      lsame = inta .eq. intb

      end
