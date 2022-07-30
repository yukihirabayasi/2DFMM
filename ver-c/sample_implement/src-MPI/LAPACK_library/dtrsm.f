      subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

      double precision alpha
      integer lda,ldb,m,n
      character diag,side,transa,uplo

      double precision a(lda,*),b(ldb,*)

      logical lsame
      external lsame

      external xerbla

      intrinsic max

      double precision temp
      integer i,info,j,k,nrowa
      logical lside,nounit,upper

      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)

      lside = lsame(side,'l')
      if (lside) then
          nrowa = m
      else
          nrowa = n
      end if
      nounit = lsame(diag,'n')
      upper = lsame(uplo,'u')

      info = 0
      if ((.not.lside) .and. (.not.lsame(side,'r'))) then
          info = 1
      else if ((.not.upper) .and. (.not.lsame(uplo,'l'))) then
          info = 2
      else if ((.not.lsame(transa,'n')) .and.  &
     &         (.not.lsame(transa,'t')) .and.  &
     &         (.not.lsame(transa,'c'))) then
          info = 3
      else if ((.not.lsame(diag,'u')) .and. (.not.lsame(diag,'n'))) then
          info = 4
      else if (m.lt.0) then
          info = 5
      else if (n.lt.0) then
          info = 6
      else if (lda.lt.max(1,nrowa)) then
          info = 9
      else if (ldb.lt.max(1,m)) then
          info = 11
      end if
      if (info.ne.0) then
          call xerbla('dtrsm ',info)
          return
      end if

      if (m.eq.0 .or. n.eq.0) return

      if (alpha.eq.zero) then
          do 20 j = 1,n
              do 10 i = 1,m
                  b(i,j) = zero
   10         continue
   20     continue
          return
      end if

      if (lside) then
          if (lsame(transa,'n')) then

              if (upper) then
                  do 60 j = 1,n
                      if (alpha.ne.one) then
                          do 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     continue
                      end if
                      do 50 k = m,1,-1
                          if (b(k,j).ne.zero) then
                              if (nounit) b(k,j) = b(k,j)/a(k,k)
                              do 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         continue
                          end if
   50                 continue
   60             continue
              else
                  do 100 j = 1,n
                      if (alpha.ne.one) then
                          do 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     continue
                      end if
                      do 90 k = 1,m
                          if (b(k,j).ne.zero) then
                              if (nounit) b(k,j) = b(k,j)/a(k,k)
                              do 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         continue
                          end if
   90                 continue
  100             continue
              end if
          else

              if (upper) then
                  do 130 j = 1,n
                      do 120 i = 1,m
                          temp = alpha*b(i,j)
                          do 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     continue
                          if (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  120                 continue
  130             continue
              else
                  do 160 j = 1,n
                      do 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          do 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     continue
                          if (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  150                 continue
  160             continue
              end if
          end if
      else
          if (lsame(transa,'n')) then

              if (upper) then
                  do 210 j = 1,n
                      if (alpha.ne.one) then
                          do 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     continue
                      end if
                      do 190 k = 1,j - 1
                          if (a(k,j).ne.zero) then
                              do 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         continue
                          end if
  190                 continue
                      if (nounit) then
                          temp = one/a(j,j)
                          do 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     continue
                      end if
  210             continue
              else
                  do 260 j = n,1,-1
                      if (alpha.ne.one) then
                          do 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     continue
                      end if
                      do 240 k = j + 1,n
                          if (a(k,j).ne.zero) then
                              do 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         continue
                          end if
  240                 continue
                      if (nounit) then
                          temp = one/a(j,j)
                          do 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     continue
                      end if
  260             continue
              end if
          else

              if (upper) then
                  do 310 k = n,1,-1
                      if (nounit) then
                          temp = one/a(k,k)
                          do 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     continue
                      end if
                      do 290 j = 1,k - 1
                          if (a(j,k).ne.zero) then
                              temp = a(j,k)
                              do 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         continue
                          end if
  290                 continue
                      if (alpha.ne.one) then
                          do 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     continue
                      end if
  310             continue
              else
                  do 360 k = 1,n
                      if (nounit) then
                          temp = one/a(k,k)
                          do 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     continue
                      end if
                      do 340 j = k + 1,n
                          if (a(j,k).ne.zero) then
                              temp = a(j,k)
                              do 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         continue
                          end if
  340                 continue
                      if (alpha.ne.one) then
                          do 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     continue
                      end if
  360             continue
              end if
          end if
      end if

      return

      end
