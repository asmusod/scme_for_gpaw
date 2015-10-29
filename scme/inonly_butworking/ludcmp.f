      subroutine ludcmp(a,n,np,indx,d)

      implicit real*8 (a-h, o-z)
C Changed by Fer
C     parameter (nmax=100,tiny=1.0e-20)
      parameter (nmax=3,tiny=1.0e-20)
      real*8 a(np,np),vv(nmax)
      integer indx(n)

C Debug
C     write(10,FMT='(1X,A)') ' Entering ludcmp'

C Debug
C     write(10,FMT='(1X,A)') ' a '
C     write(10,FMT='(1X,3E23.15)') (a(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (a(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (a(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' n, np, indx '
C     write(10,FMT='(1X,5I5)') n, np, indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' d '
C     write(10,FMT='(1X,F20.14)') d
C     write(10,FMT='(1X,A)') ' '

      d=1.d0
      do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
            if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11      continue
         if (aamax.eq.0.) then 
            print '(A)', 'singular matrix in LUDCMP. Cannot shake.'
            stop
         end if
         vv(i)=1.d0/aamax
 12   continue

C Debug
C     write(10,FMT='(1X,A)') ' vv '
C     write(10,FMT='(1X,3E23.15)') vv
C     write(10,FMT='(1X,A)') ' '

      do 19 j=1,n
         if (j.gt.1) then
            do 14 i=1,j-1
               sum=a(i,j)
               if (i.gt.1)then
                  do 13 k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
 13               continue
                  a(i,j)=sum
               endif
 14         continue
         endif
         aamax=0.d0
         do 16 i=j,n
            sum=a(i,j)
            if (j.gt.1)then
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
 15            continue
               a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
 16      continue
         if (j.ne.imax)then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
 17         continue
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(j.ne.n)then
            if(a(j,j).eq.0.)a(j,j)=tiny
            dum=1.d0/a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
 18         continue
         endif
 19   continue
      if(a(n,n).eq.0.)a(n,n)=tiny

C Debug
C     write(10,FMT='(1X,A)') ' a '
C     write(10,FMT='(1X,3E23.15)') (a(1,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (a(2,p),p=1,3)
C     write(10,FMT='(1X,3E23.15)') (a(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' n, np, indx '
C     write(10,FMT='(1X,5I5)') n, np, indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' d '
C     write(10,FMT='(1X,F20.14)') d
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' Leaving ludcmp'

      return
      end
c-------------------------------------------------------------------
      subroutine lubksb(a,n,np,indx,b)
      real*8 a(np,np),b(n)
      integer indx(n)

C Debug
C     write(10,FMT='(1X,A)') ' Entering lubksb --------------------'

C Debug
C     write(10,FMT='(1X,A)') ' n np '
C     write(10,FMT='(1X,2I5)') n, np
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' a '
C     write(10,FMT='(1X,3F20.14)') (a(1,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (a(2,p),p=1,3)
C     write(10,FMT='(1X,3F20.14)') (a(3,p),p=1,3)
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' indx '
C     write(10,FMT='(1X,3I5)') indx
C     write(10,FMT='(1X,A)') ' '

C Debug
C     write(10,FMT='(1X,A)') ' b '
C     write(10,FMT='(1X,3F20.14)') b
C     write(10,FMT='(1X,A)') ' '

      ii=0
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 11         continue
         else if (sum.ne.0.) then
            ii=i
         endif
         b(i)=sum
 12   continue
      do 14 i=n,1,-1
         sum=b(i)
         if(i.lt.n)then
            do 13 j=i+1,n
               sum=sum-a(i,j)*b(j)
 13         continue
         endif
         b(i)=sum/a(i,i)
 14   continue

C Debug
C     write(10,FMT='(1X,A)') ' Leaving lubksb --------------------'

      return
      end
