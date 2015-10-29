      subroutine tangentialVel(r1, r2, r3, v1, v2, v3, w, mMolecule)

      implicit real*8 (a-h,o-z)
      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'

      real*8 r1(3), r2(3), r3(3), v1(3), v2(3), v3(3), rCM(3)
      real*8 x1(3), x2(3), x3(3), w(3), xx1, xx2, mMolecule, wCrystal(3)

c     Recover the molecule in case that it was broken due to the
c     periodic boundary conditions
      call reconstruct(r1, r2, r3)

c     Find the unit vectors along the inertia principal axis of the
c     molecule, x1, x2, x3 
      do i = 1, 3
         x1(i) = r1(i) - r2(i)
         x2(i) = (r3(i)-r1(i)) + (r3(i)-r2(i))
      end do

      xx1 = sqrt(dot(x1,x1))
      xx2 = sqrt(dot(x2,x2))
      do i = 1, 3
         x1(i) = x1(i) / xx1
         x2(i) = x2(i) / xx2
      end do

      call cross(x1, x2, x3)

c     Find the center of mass, rCM
      do i = 1, 3
         rCM(i) = (amass(1)* r1(i) + amass(1)* r2(i) + 
     $        amass(2)* r3(i))/ mMolecule
      end do

c     Change the coordinates to a frame at the center of mass
      do i = 1, 3
         r1(i) = r1(i) - rCM(i)
         r2(i) = r2(i) - rCM(i)
         r3(i) = r3(i) - rCM(i)
      end do

c     Transform the angular velocity w to the reference frame of the crystal.
      do i = 1, 3
         wCrystal(i) = x1(i) * w(1) + x2(i) * w(2) + x3(i) * w(3)
      end do

c     Compute tangential velocities, r x w
      call cross(r1, wCrystal, v1)
      call cross(r2, wCrystal, v2)
      call cross(r3, wCrystal, v3)

      return
      end

c-----------------------------------------------------------------------
      subroutine reconstruct(r1, r2, r3)

c     This routine recovers finds the minimum image of each hydrogen in
c     a molecule, given the oxigen.

      implicit real*8 (a-h,o-z)
      include '../commonblks/comgeom.cmn'    

      real*8 r1(3), r2(3), r3(3), d, s

c     Check minimum image in x direction
      d = r3(1) - r1(1)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. 1.d0/2.d0) then
         r1(1) = r1(1) + s
      end if
      d = r3(1) - r2(1)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. 1.d0/2.d0) then
         r2(1) = r2(1) + s
      end if

c     Check minimum image in y direction
      d = r3(2) - r1(2)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. ay/ax/2.d0) then
         r1(2) = r1(2) + s * ay/ax
      end if
      d = r3(2) - r2(2)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. ay/ax/2.d0) then
         r2(2) = r2(2) + s * ay/ax
      end if

c     Check minimum image in z direction
      d = r3(3) - r1(3)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. az/ax/2.d0) then
         r1(3) = r1(3) + s * az/ax
      end if
      d = r3(3) - r2(3)
      s = SIGN(1.0D+00,d)
      if (abs(d) .gt. az/ax/2.d0) then
         r2(3) = r2(3) + s * az/ax
      end if

      return
      end





