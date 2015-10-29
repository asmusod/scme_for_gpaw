c----------------------------------------------------------------------+
c     Routine to fix molecules broken due to the periodic boundary     |
c     conditions                                                       |
c----------------------------------------------------------------------+
      subroutine recoverMolecules(raOri, ra, nH, nO, a, a2)

      implicit real*8 (a-h,o-z)
      include '../commonblks/parameters.cmn'

      real*8 raOri(maxCoo), a(3), a2(3), ra(maxCoo), dist
      integer nH,nO

      do i = 1, nO
         do l = 1, 3
            do n = 1, 2
               index = l+3*(n-1)+6*(i-1)
               dist = raOri(index) - raOri(l+3*(i-1+nH))
               if (dist .gt. a2(l)) then
                  ra(index) = raOri(index) - a(l)
               elseif(dist .lt. -a2(l)) then
                  ra(index) = raOri(index) + a(l)
               else
                  ra(index) = raOri(index)
               end if
            end do
            ra(l+3*(i-1+nH)) = raOri(l+3*(i-1+nH))
         end do
      end do
      return
      end
c----------------------------------------------------------------------+
c     Routine to calculate the center of mass of each molecule         |
c----------------------------------------------------------------------+
      subroutine calcCentersOfMass(ra, nM, rCM)

      implicit real*8 (a-h,o-z)
      include '../commonblks/parameters.cmn'
      integer nM, iH1, iH2, iO, i, j
      real*8 rCM(3,maxCoo/3), ra(maxCoo)

      do i = 1, nM
         iH2 = 2 * i
         iH1 = iH2 - 1
         iO  = 2 * nM + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO  = 3 * (iO - 1)
         do j = 1, 3
            rCM(j,i) = (ra(iH1+j) + ra(iH2+j) + 16.d0 * ra(iO+j)) / 
     $           18.d0
         end do
      end do

      return
      end
c----------------------------------------------------------------------+
c     Routine to calculate the principal axes of each molecule         |
c----------------------------------------------------------------------+
      subroutine findPpalAxesOLD(ra, nM, x)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      integer nM, iH1, iH2, iO, i, j
      real*8 x(3,3,maxCoo/3), ra(maxCoo), r11, r21

      do i = 1, nM
         iH2 = 2 * i
         iH1 = iH2 - 1
         iO  = 2 * nM + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO  = 3 * (iO - 1)
         do j = 1, 3
            x(j,1,i) = -(ra(iH1+j) + ra(iH2+j) - 2.d0 * ra(iO+j))
            x(j,2,i) = ra(iH2+j) - ra(iH1+j)
         end do
         r11 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i)*
     $        x(3,1,i))
         r21 = sqrt(x(1,2,i)*x(1,2,i) + x(2,2,i)*x(2,2,i) + x(3,2,i)*
     $        x(3,2,i))
         do j = 1, 3
            x(j,1,i) = x(j,1,i) / r11
            x(j,2,i) = x(j,2,i) / r21
         end do
         x(1,3,i) = x(2,1,i) * x(3,2,i) - x(3,1,i) * x(2,2,i)
         x(2,3,i) = x(3,1,i) * x(1,2,i) - x(1,1,i) * x(3,2,i)
         x(3,3,i) = x(1,1,i) * x(2,2,i) - x(2,1,i) * x(1,2,i)
      end do
      return
      end
c----------------------------------------------------------------------+
c     Routine to calculate the principal axes of each molecule         |
c----------------------------------------------------------------------+
      subroutine findPpalAxes(ra, nM, x)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      integer nM, iH1, iH2, iO, i, j
      real*8 x(3,3,maxCoo/3), ra(maxCoo), r11, r21

C Debug
      integer*4 p

      do i = 1, nM
         iH2 = 2 * i
         iH1 = iH2 - 1
         iO  = 2 * nM + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO  = 3 * (iO - 1)
         do j = 1, 3
            x(j,3,i) = -(ra(iH1+j) + ra(iH2+j) - 2.d0 * ra(iO+j))
            x(j,1,i) = ra(iH2+j) - ra(iH1+j)
         end do
         r11 = sqrt(x(1,3,i)*x(1,3,i) + x(2,3,i)*x(2,3,i) + x(3,3,i)*
     $        x(3,3,i))
         r21 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i)*
     $        x(3,1,i))

         do j = 1, 3
            x(j,3,i) = x(j,3,i) / r11
            x(j,1,i) = x(j,1,i) / r21
         end do
         x(1,2,i) = x(2,3,i) * x(3,1,i) - x(3,3,i) * x(2,1,i)
         x(2,2,i) = x(3,3,i) * x(1,1,i) - x(1,3,i) * x(3,1,i)
         x(3,2,i) = x(1,3,i) * x(2,1,i) - x(2,3,i) * x(1,1,i)
      end do

C Debug
C Print the principal axes matrix for each molecule
C     do i=1,nM
C       write(10,FMT='(A,I4)') ' Axes for Molecule: ', i
C       write(10,FMT='(3F10.5)') (x(1,p,i),p=1,3)
C       write(10,FMT='(3F10.5)') (x(2,p,i),p=1,3)
C       write(10,FMT='(3F10.5)') (x(3,p,i),p=1,3)
C     end do

      return
      end

c----------------------------------------------------------------------+
c     Routine to rotate the multipoles to the orientation of the       |
c     molecule                                                         |
c----------------------------------------------------------------------+
      subroutine rotatePoles(d0, q0, o0, h0, dpole, qpole, opole, hpole,
     $     nM, x) 

      implicit none

      include '../commonblks/parameters.cmn'
      integer nM, i, j, k, l, ii, jj, kk, ll, m
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
      
      real*8 x(3,3,maxCoo/3)
      real*8 d0(3), q0(3,3), o0(3,3,3), h0(3,3,3,3)
      
      do m = 1, nM
         do l = 1, 3
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     hpole(i,j,k,l,m) = 0.d0
                  end do
                  opole(j,k,l,m) = 0.d0
               end do
               qpole(k,l,m) = 0.d0
            end do
            dpole(l,m) = 0.d0
         end do
      end do
      
      do m = 1, nM
         do l = 1, 3
            do ll = 1, 3
               do k = 1, 3
                  do kk = 1, 3
                     do j = 1, 3
                        do jj = 1, 3
                           do i = 1, 3
                              do ii = 1, 3
                                 hpole(i,j,k,l,m) = hpole(i,j,k,l,m) +
     $                                x(i,ii,m) * x(j,jj,m) * x(k,kk,m)
     $                                * x(l,ll,m) * h0(ii,jj,kk,ll) 
                              end do
                           end do
                           opole(j,k,l,m) = opole(j,k,l,m) + x(j,jj,m)
     $                          * x(k,kk,m)* x(l,ll,m) * o0(jj,kk,ll) 
                        end do
                     end do
                     qpole(k,l,m) = qpole(k,l,m) + x(k,kk,m) * x(l,ll,m)
     $                    * q0(kk,ll)
                  end do
               end do
               dpole(l,m) = dpole(l,m) + x(l,ll,m) * d0(ll)
            end do
         end do
      end do

C     do i=1,3
C       write(10,FMT='(1X,3F10.5)') (x(i,j,1),j=1,3)
C     end do

C     do i=1,3
C       write(10,FMT='(1X,3F10.5)') (x(i,j,2),j=1,3)
C     end do
      
      return
      end
c----------------------------------------------------------------------+
c     Routine to rotate the multipoles to the orientation of the       |
c     molecule                                                         |
c----------------------------------------------------------------------+
      subroutine setUnpolPoles(dpole, qpole, dpole0, qpole0, nM) 

      implicit none

      include '../commonblks/parameters.cmn'
      integer nM, i, j, k
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)

      do i = 1, nM
         do j = 1, 3
            dPole(j,i) = dpole0(j,i)
            do k = 1, 3
               qpole(k,j,i) = qpole0(k,j,i)
            end do
         end do
      end do
      return
      end
c----------------------------------------------------------------------+
c     Routine to rotate the multipoles to the orientation of the       |
c     molecule                                                         |
c----------------------------------------------------------------------+
      subroutine rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, 
     $     x)

      implicit none

      include '../commonblks/parameters.cmn'
      integer nM, i, j, k, l, ii, jj, kk, ll, m

      real*8 x(3,3,maxCoo/3)
      real*8 hp(3,3,3,maxCoo/3), dq(3,3,3,maxCoo/3), dd(3,3,maxCoo/3)
      real*8 qq(3,3,3,3,maxCoo/3)
      real*8 hp0(3,3,3), dq0(3,3,3), dd0(3,3), qq0(3,3,3,3)

      do m = 1, nM
         do k = 1, 3
            do j = 1, 3
               do i = 1, 3
                  do l = 1, 3
                     qq(l,i,j,k,m) = 0.d0
                  end do
                  dq(i,j,k,m) = 0.d0
                  hp(i,j,k,m) = 0.d0
               end do
               dd(j,k,m) = 0.d0
            end do
         end do
      end do
      
      do m = 1, nM
         do k = 1, 3
            do kk = 1, 3
               do j = 1, 3
                  do jj = 1, 3
                     do i = 1, 3
                        do ii = 1, 3
                           do l = 1, 3
                              do ll = 1, 3
                                 qq(l,i,j,k,m) = qq(l,i,j,k,m) + x(l,ll
     $                                ,m) *  x(i,ii,m) * x(j,jj,m)* x(k
     $                                ,kk,m) * qq0(ll,ii,jj,kk) 

                              end do
                           end do
                           dq(i,j,k,m) = dq(i,j,k,m) + x(i,ii,m) * x(j
     $                          ,jj,m)* x(k,kk,m) * dq0(ii,jj,kk) 
                           hp(i,j,k,m) = hp(i,j,k,m) + x(i,ii,m) * x(j
     $                          ,jj,m)* x(k,kk,m) * hp0(ii,jj,kk)
                        end do
                     end do
                     dd(j,k,m) = dd(j,k,m) + x(j,jj,m) * x(k,kk,m)
     $                    * dd0(jj,kk) 
                  end do
               end do
            end do
         end do
      end do

      return
      end
c----------------------------------------------------------------------+
c     Add all the fields                                               |
c----------------------------------------------------------------------+
      subroutine addFields(eH, eD, eT, nM)

      implicit none
      include '../commonblks/parameters.cmn'
      integer i, j, nM
      real*8 eH(3,maxCoo/3), eD(3,maxCoo/3) 
      real*8 eT(3,maxCoo/3)

      do i = 1, nM
         do j = 1, 3
            eT(j,i) = eH(j,i) + eD(j,i)
         end do
      end do
      return
      end
c----------------------------------------------------------------------+
c     Add the derivative of all the fields                             |
c----------------------------------------------------------------------+
      subroutine addDfields(dEhdr, dEddr, dEtdr, nM)

      implicit none
      include '../commonblks/parameters.cmn'
      integer i, j, k, nM
      real*8 dEhdr(3,3,maxCoo/3)
      real*8 dEtdr(3,3,maxCoo/3), dEddr(3,3,maxCoo/3)

      do i = 1, nM
         do j = 1, 3
            do k = 1, 3
               dEtdr(k,j,i) = dEhdr(k,j,i) + dEddr(k,j,i)
c               dEtdr(j,k,i) = dEhdr(k,j,i)
            end do
         end do
         do j = 2, 3
            do k = 1, j-1
               dEtdr(j,k,i) = dEtdr(k,j,i)
            end do
         end do
      end do
      return
      end
c-----------------------------------------------------------------------      
      subroutine applyPBC(r, sep, j)

      implicit real*8 (a-h,o-z)
      include '../commonblks/comgeom.cmn'
      real*8 r(3), a(3), sep(3)
      integer j

c     Size of the simulation cell
      a(1) = ax/2.0d0
      a(2) = ay/2.0d0
      a(3) = az/2.0d0

      do i = 1, 3
c         if (j*sep(i) .gt. a(i)) then
         if (r(i) .gt. a(i)) then
            r(i) = r(i) - 2.0d0 * a(i)
c         elseif (j*sep(i) .lt. -a(i)) then
         elseif (r(i) .lt. -a(i)) then
            r(i) = r(i) + 2.0d0 * a(i)
         end if
      end do

      j = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine SF(r, swFunc)

      implicit none
      real*8 r, swFunc, x, x2, x3, rSW, rCut
      real*8 rL1, rL2, rH1, rH2, dr

c      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
      data rL1, rH1, rL2, rH2 / 0.d0, 5.d0, 9.d0, 11.d0 /
c      data rL1, rH1, rL2, rH2 / 0.d0, 5.d0, 11.d0, 13.d0 /
      save 


c     Kroes
c$$$      rSW = 8.d0
c$$$      rCut = 9.d0
c$$$
c$$$      x = (r - rSW)/(rCut - rSW)
c$$$
c$$$      if (r. lt. rSW) then
c$$$         swFunc = 1.0d0
c$$$      else if(r .lt. rCut) then
c$$$         x2 = x * x
c$$$         x3 = x2 * x
c$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
c$$$      else
c$$$         swFunc = 0.0d0
c$$$      end if

      if ((r .ge. rH2) .or. (r .le. rL1)) then
         swFunc = 0.0d0
      else if ((r .ge. rH1) .and. (r .le. rL2)) then
         swFunc = 1.0d0
      else if (r .lt. rH1) then

         call tang_toenniesN(r, swFunc, 6)

c$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
c$$$         x2 = x * x
c$$$         x3 = x2 * x
c$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
      else
         x = (r - rL2)/(rH2 - rL2)
         x2 = x * x
         x3 = x2 * x
         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
      end if

c      swFunc = 1.d0

      return
      end

c-----------------------------------------------------------------------
      subroutine SFdsf(r, swFunc, dSdr)

      implicit none
      real*8 r, swFunc, x, x2, x3, dSdr, rSW, rCut
      real*8 rL1, rL2, rH1, rH2, dr

c      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
      data rL1, rH1, rL2, rH2 / 0.d0, 5.d0, 9.d0, 11.d0 /
c      data rL1, rH1, rL2, rH2 / 0.d0, 5.d0, 11.d0, 13.d0 /
      save 

c     Kroes
c$$$      rSW = 8.d0
c$$$      rCut = 9.d0


c$$$      x = (r - rSW)/(rCut - rSW)
c$$$
c$$$      if (r. lt. rSW) then
c$$$         swFunc = 1.0d0
c$$$         dSdr = 0.0d0
c$$$      else if(r .lt. rCut) then
c$$$         x2 = x * x
c$$$         x3 = x2 * x
c$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
c$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) / (rCut-rSW)
c$$$      else
c$$$         swFunc = 0.0d0
c$$$         dSdr = 0.0d0
c$$$      end if


      if ((r .ge. rH2) .or. (r .le. rL1)) then
         swFunc = 0.0d0
         dSdr = 0.0d0
      else if ((r .ge. rH1) .and. (r .le. rL2)) then
         swFunc = 1.0d0
         dSdr = 0.0d0
      else if (r .lt. rH1) then

cc         call switchCloseDsf(r, swFunc, dSdr)
         call tang_toenniesNdF(r, swFunc, dSdr, 6)

c$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
c$$$         dr = - 1.d0 / (rH1 - rL1)
c$$$         x2 = x * x
c$$$         x3 = x2 * x
c$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
c$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) * dr
      else
         x = (r - rL2)/(rH2 - rL2)
         dr = 1.d0 / (rH2 - rL2)
         x2 = x * x
         x3 = x2 * x
         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) * dr
      end if

c      swFunc = 1.d0
c      dSdr = 0.d0

      return
      end

c-----------------------------------------------------------------------      
c----------------------------- The End ---------------------------------
c-----------------------------------------------------------------------      
