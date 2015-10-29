      subroutine atomicForces(fCM, tau, ra, rCM, nM, fa)

      implicit none
      include '../commonblks/parameters.cmn'

      integer nM
      real*8 fCM(3,maxCoo/3), rCM(3,maxCoo/3), tau(3,maxCoo/3)
      real*8 ra(maxCoo), fa(maxCoo)

      integer n, iH1, iH2, iO, i
      real*8 fTot(3), tt(3), r1(3), r2(3), r3(3), f1(3), f2(3), f3(3) 
      integer flag

      do n = 1, nM
         iH1 = 6*(n-1)
         iH2 = iH1 + 3
         iO = 3 * (2*nM + n - 1)

         do i = 1, 3
            fTot(i) = fCM(i,n)
            tt(i) = tau(i,n)
            r1(i) = ra(iH1+i) - rCM(i,n)
            r2(i) = ra(iH2+i) - rCM(i,n)
            r3(i) = ra(iO+i)  - rCM(i,n)
         end do
c     as comment out, it's only correcting the force 
c     with regards to RATTLE constraints:
         call molForce3(f1, f2, f3, r1, r2, r3, fTot, tt, flag)
c$$$         if (flag .gt. 0) then
c$$$            print *, n
c$$$         end if

         do i = 1, 3
            fa(iH1+i) = f1(i)
            fa(iH2+i) = f2(i)
            fa(iO+i)  = f3(i)
         end do
      end do

      return
      end
