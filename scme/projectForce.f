      subroutine projectForce
      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/constraints.cmn'

      real*8 c, cDf

      c = 0.d0
      cDf = 0.d0
c      cdv = 0.d0
      do i = 1, 3*natoms
         cDf = cDf + constraint(i)*fa(i)
c         cDv = cDv + constraint(i)*va(i)
      end do

      do i = 1, 3*natoms
         fa(i) = fa(i) - cDf * constraint(i)
c         va(i) = va(i) - cDv * constraint(i)
      end do

      return
      end

c-----------------------------------------------------------------------
c     General constraint of the centers of mass

      subroutine projectForce2
      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/constraints.cmn'

      real*8 fcm(maxatoms*3)
c
c     Project forces
c     
      do i = 1, nAtms(2)
         iH1 = 2 * i - 1
         iH2 = iH1 + 1
         iO = nAtms(1) + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO = 3 * (iO - 1)
         ic = 3 * (i-1)
         do j = 1, 3
            fcm(ic+j) = fa(iH1+j) + fa(iH2+j) + fa(iO+j)
         end do
      end do

      cDf = 0.d0
      do i = 1, 3 * nAtms(2)         
         cDf = cDf + fcm(i) * constraint(i)
      end do

      do i = 1, 3 * nAtms(2)         
         fcm(i) = cDf * constraint(i)
      end do

      do i = 1, nAtms(2)
         iH1 = 2 * i - 1
         iH2 = iH1 + 1
         iO = nAtms(1) + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO = 3 * (iO - 1)
         ic = 3 * (i-1)
         do j = 1, 3
            fa(iH1+j) = fa(iH1+j) - fcm(ic+j) / 18.
            fa(iH2+j) = fa(iH2+j) - fcm(ic+j) / 18.
            fa(iO+j)  = fa(iO+j)  - fcm(ic+j) * 16. / 18.
         end do
      end do

c
c     Project velocities
c     
      do i = 1, nAtms(2)
         iH1 = 2 * i - 1
         iH2 = iH1 + 1
         iO = nAtms(1) + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO = 3 * (iO - 1)
         ic = 3 * (i-1)
         do j = 1, 3
            fcm(ic+j) = va(iH1+j) + va(iH2+j) + 16.d0 * va(iO+j)
            fcm(ic+j) = fcm(ic+j) / 18.d0
         end do
      end do

      cDv = 0.d0
      do i = 1, 3 * nAtms(2)         
         cDv = cDv + fcm(i) * constraint(i)
      end do

      do i = 1, 3 * nAtms(2)         
         fcm(i) = cDv * constraint(i)
      end do

      do i = 1, nAtms(2)
         iH1 = 2 * i - 1
         iH2 = iH1 + 1
         iO = nAtms(1) + i

         iH1 = 3 * (iH1 - 1)
         iH2 = 3 * (iH2 - 1)
         iO = 3 * (iO - 1)
         ic = 3 * (i-1)
         do j = 1, 3
            va(iH1+j) = va(iH1+j) - fcm(ic+j)
            va(iH2+j) = va(iH2+j) - fcm(ic+j)
            va(iO+j)  = va(iO+j)  - fcm(ic+j)
         end do
      end do

      return
      end
