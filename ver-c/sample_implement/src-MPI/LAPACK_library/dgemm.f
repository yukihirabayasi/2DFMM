      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

      double precision alpha,beta
      integer k,lda,ldb,ldc,m,n
      character transa,transb

      double precision a(lda,*),b(ldb,*),c(ldc,*)

      logical lsame
      external lsame

      external xerbla

      intrinsic max

      double precision temp
      integer i,info,j,l,ncola,nrowa,nrowb
      logical nota,notb

      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)

      nota = lsame(transa,'n')
      notb = lsame(transb,'n')
      if (nota) then
          nrowa = m
          ncola = k
      else
          nrowa = k
          ncola = m
      end if
      if (notb) then
          nrowb = k
      else
          nrowb = n
      end if

      info = 0
      if ((.not.nota) .and. (.not.lsame(transa,'c')) .and.  &
   &      (.not.lsame(transa,'t'))) then
          info = 1
      else if ((.not.notb) .and. (.not.lsame(transb,'c')) .and.  &
   &           (.not.lsame(transb,'t'))) then
          info = 2
      else if (m.lt.0) then
          info = 3
      else if (n.lt.0) then
          info = 4
      else if (k.lt.0) then
          info = 5
      else if (lda.lt.max(1,nrowa)) then
          info = 8
      else if (ldb.lt.max(1,nrowb)) then
          info = 10
      else if (ldc.lt.max(1,m)) then
          info = 13
      end if
      if (info.ne.0) then
          call xerbla('dgemm ',info)
          return
      end if

      if ((m.eq.0) .or. (n.eq.0) .or.  &
   &      (((alpha.eq.zero).or. (k.eq.0)).and. (beta.eq.one))) return

      if (alpha.eq.zero) then
          if (beta.eq.zero) then
              do 20 j = 1,n
                  do 10 i = 1,m
                      c(i,j) = zero
   10             continue
   20         continue
          else
              do 40 j = 1,n
                  do 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             continue
   40         continue
          end if
          return
      end if

      if (notb) then
          if (nota) then

              do 90 j = 1,n
                  if (beta.eq.zero) then
                      do 50 i = 1,m
                          c(i,j) = zero
   50                 continue
                  else if (beta.ne.one) then
                      do 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 continue
                  end if
                  do 80 l = 1,k
                      temp = alpha*b(l,j)
                      do 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 continue
   80             continue
   90         continue
          else

              do 120 j = 1,n
                  do 110 i = 1,m
                      temp = zero
                      do 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 continue
                      if (beta.eq.zero) then
                          c(i,j) = alpha*temp
                      else
                          c(i,j) = alpha*temp + beta*c(i,j)
                      end if
  110             continue
  120         continue
          end if
      else
          if (nota) then

              do 170 j = 1,n
                  if (beta.eq.zero) then
                      do 130 i = 1,m
                          c(i,j) = zero
  130                 continue
                  else if (beta.ne.one) then
                      do 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 continue
                  end if
                  do 160 l = 1,k
                      temp = alpha*b(j,l)
                      do 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  150                 continue
  160             continue
  170         continue
          else

              do 200 j = 1,n
                  do 190 i = 1,m
                      temp = zero
                      do 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 continue
                      if (beta.eq.zero) then
                          c(i,j) = alpha*temp
                      else
                          c(i,j) = alpha*temp + beta*c(i,j)
                      end if
  190             continue
  200         continue
          end if
      end if

      return

      end
