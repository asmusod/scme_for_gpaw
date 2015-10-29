c     New version of nudgedEB for all the molecules in chain. This
c     routine has been optimized from the previous one nudgeEB.f 
c     making the memory requirement 38 times smaller. 

      subroutine SADDLEFIND

c     nudgedEB.f:
c     in version f the dwitching function goes from 0 to 1 in the 
c     range between cosPhilarge and cosPhismal  
c     perameters comng from the combaths commonblk.
c     e: Implement switching fucntion to gradually turn on
c     perpendicular components
c     of the spring force when paths get kinky:

c       f(x) = 0                             if cosphi>cosphismall
c            = 1                             if cosphi<cosphilarge
c             switch = 0.5*(1.0+cos(Pi*(djdotdjn-cosPhilarge)/
c                                   (cosPhismall-cosPhilarge)))   else
c          where cosPhimax which must be > 0

c     This routine nudges the elastic band by modifying the forces.  
c     The elastic band is represented as an FPI
c     fixed-endpoint chain with the ends hopefully on the minimum
c     -energy path. The path vector is ds, its length is unimportant
c     and its direction is at half the bending angle between the
c     adjacent images.  Then, the force due to the potential and the
c     force due to the springs are projected into components 
c     parallel and perendicular to ds.  
c     The perpendicular component of the real force, and the
c     parallel component of the spring force, are kept and added
c     together to find the new effective force.
c     This routine is called in place of FPIforce.

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comoutcntrl.cmn'
      include '../commonblks/combaths.cmn'

      data Pi / 3.14159265358979312d0 /
      real*8 fs(maxcoo), ds(maxcoo), dj(maxcoo), djn(maxcoo), switch
      real*8 box(3)

c     Now, ds is found by first finding the bending angle phi between
c     the previous, current, and next images; i.e. the vectors r1 and
c     r2.  ds is bent from r1 by an angle of phi/2 in the direction of
c     r2. Then take the projections of the forces due to U and the
c     springs onto ds.

c     There is a maximum bend per image (Phi) allowed for the path.
c     Phi at image j is Pi, minus the angle made by the points Rj-1, Rj
c     , and Rj+1.  The cosine of Phi is the dot product of the
c     normalized vectors dj and djn.  If the bend angle is greater than
c     a cutoff Phimax, the perpendicular component of the
c     spring force is not thrown out, which allows the springs to
c     straighten the path.

      box(1) = ax
      box(2) = ay
      box(3) = az

      nat = natoms/nimFPI(1)
      do j = 2, nimFPI(1)-1
         l = 0
         fuperpTot2 = 0.d0
         fa2 = 0.d0
         fs2 = 0.d0

         do i = 1, nat
            index = 3 * nimFPI(1) * (i-1) + 3 * (j-1)
            do k = 1, 3
               l = l + 1
c     dj: displacement vector to previous image
               dj(l) = ax * (ra(index-3+k) - ra(index+k))
c     djn: displacement vector to next image
               djn(l) = ax * (ra(index+3+k) - ra(index+k))
c     apply periodic boundary conditions to dj and djn
               if (dj(l) .lt. -box(k)/2.) then
                  dj(l) = dj(l) + box(k)
               else if (dj(l) .gt. box(k)/2.) then
                  dj(l) = dj(l) - box(k)
               end if
               if (djn(l) .lt. -box(k)/2.) then
                  djn(l) = djn(l) + box(k)
               else if (djn(l) .gt. box(k)/2.) then
                  djn(l) = djn(l) - box(k)
               end if
c     calculate spring force
               fs(l) = sprconFPI(1) * (dj(l) + djn(l)) !spring force
c     find ds: vector from previous to next image
               ds(l) = djn(l) - dj(l)
            end do
         end do

c     calculate the length of ds, dj and djn
         dslen = 0.
         dj2   = 0.
         djn2  = 0.
         do i = 1, 3*nat
            dslen = dslen + ds(i)*ds(i)
            dj2   = dj2 + dj(i)*dj(i)
            djn2  = djn2 + djn(i)*djn(i)
         end do
         dslen = sqrt(dslen)
         dj1   = sqrt(dj2)
         djn1  = sqrt(djn2)

c     Normalize ds, dj & djn, and calculate the components of the spring
c     and potnetial forces along the path.
         djdotdjn = 0.
         fudotds  = 0.
         fsdotds  = 0.
         l = 0
         do i = 1, nat
            index = 3 * nimFPI(1) * (i-1) + 3 * (j-1)
            do k = 1, 3
               l = l + 1
               dj(l) = dj(l)/dj1
               djn(l) = djn(l)/djn1
               ds(l) = ds(l)/dslen

               djdotdjn = djdotdjn + dj(l)*djn(l)
               fudotds = fudotds + FA(index+k)*ds(l)
               fsdotds = fsdotds + fs(l)*ds(l)
            end do
         end do

         cosPhi = -djdotdjn
         switch = 0.0
         if (cosPhi .lt. cosPhismall) then
            if (cosPhi .gt. cosPhilarge) then
c     In this case the path is too kinky and some of the spring force 
c     perpendicular to the path is kept, just write out a warning here
c     and evaluate the switching function: 

               switch = 0.5 * (1.0 + cos(Pi * (cosPhi-cosPhilarge)/
     $              (cosPhismall - cosPhilarge)))
            else
               switch = 1.0
            end if
            if(Iqkmin .and. switch .gt. 0.5) 
     +           write(lunqckmin,301) im,djdotdjn,switch
 301        format('                                    ',
     +           '** spfn: Kinky path condition at image ',i4,
     +           ', djdotdjn =',g12.4,' switch =',g12.4)
         endif

c     We now have the (normalized) vector ds, and the dot products
c     fu*ds and fs*ds.  Now it's time to do some force projections.
c     We will keep the perpendicular component of fu (the gradient) and the
c     parallel component of fs (the spring force).
              
         l = 0
         do i = 1, nat
            index = 3 * nimFPI(1) * (i-1) + 3 * (j-1)
            do k = 1, 3
               l = l + 1
               fuperp = FA(index+k) - fudotds * ds(l)

c     Calculate the magnitude of the perpendicular component of the
c     total force.
               fuperpTot2 = fuperpTot2 + fuperp * fuperp
               fa2 = fa2 + FA(index+k) * fa(index+k)
cERB               fs2 = fs2 + fs(l) * fs(l)

               fspar = fsdotds * ds(l)
               fsperp = fs(l) - fspar

               FA(index+k) = fuperp + fspar + switch * fsperp
            end do
         end do
         fuperpT2(j) = fuperpTot2
         faT2(j) = fa2
      end do                    ! loop over images

      return
      end
