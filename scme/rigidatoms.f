c   Rigid atoms:
c      this subroutine zeros the velocity and force of
c      permafrost atoms, endatoms in FPI chains with tag 100-399, and the
c      appropriate components of constrained atoms (tags 1 and 4-9)
c      and images in constrained FPI chains (with tags 500-599).
c  !!!!! temperature gets messed up when there are many FPI endatoms, check
c        Gregs solution to this.           

      subroutine rigidatoms(v2perm)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comluns.cmn'
      integer chainLength
      real*8 vCM(3), fTotal(3), n(3)

      real*8 rCMi(3), rCMk(3), dr(3), ddd
c-----------------------------------------------------------------------
c Permafrost, tag=1:
                                                                      
c     Zero both the velocity and force on all permafrost atoms:
      do 136 ip=1,naperm
         ik=iperm(ip)
c$$$         v2perm=v2perm+va(3*(ik-1)+1)**2+va(3*(ik-1)+2)**2+
c$$$     +        va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
c     from permafrost atoms.
         va(3*(ik-1)+1)=0.0
         va(3*(ik-1)+2)=0.0
         va(3*(ik-1)+3)=0.0
         fa(3*(ik-1)+1)=0.0
         fa(3*(ik-1)+2)=0.0
         fa(3*(ik-1)+3)=0.0
 136  continue

cd           write(lunout,*)  '  RIGIDATOMS:  naperm = ',naperm
cd           write(lunout,*) '  After permafrost operation: time=',time
cd           write(lunout,233) 1,(va(3*k),k=1,natms(1))
cd           write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))
cd233        format('   z component of velocities for type',i1,'atoms:'/
cd   +                 5g15.6/5g15.6/5g15.6/5g15.6)

c-----------------------------------------------------------------------
c  FPI chains:
                                                                            
      if(nFPI .gt. 0) then
         do 137 ip=1,nFPI
c     loop over FPI chains:

c     Fixed-end FPI:
                                                                      
            if((itagFPI(ip) .ge. 100 .and. itagFPI(ip) .lt. 200) .or.
     +           (itagFPI(ip) .ge. 500 .and. itagFPI(ip) .lt. 600)) 
     +           then
c for tags 100-199 and 500-599 zero the velocity and force on the end
c atoms:

               ik=iFPI(ip)
c     ik points at the first atom in FPI chain ip.
                               
c$$$               v2perm=v2perm+va(3*(ik-1)+1)**2+va(3*(ik-1)+2)**2+
c$$$     +              va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
c     from permafrost atoms.
                                          
               va(3*(ik-1)+1)=0.0
               va(3*(ik-1)+2)=0.0
               va(3*(ik-1)+3)=0.0
               fa(3*(ik-1)+1)=0.0
               fa(3*(ik-1)+2)=0.0
               fa(3*(ik-1)+3)=0.0

               ik=iFPI(ip)+nimFPI(ip)-1
c     ik points at the last atom in FPI chain ip.
                
c$$$               v2perm=v2perm+va(3*(ik-1)+1)**2+va(3*(ik-1)+2)**2+
c$$$     +              va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
c     from atoms.
                                                     
               va(3*(ik-1)+1)=0.0
               va(3*(ik-1)+2)=0.0
               va(3*(ik-1)+3)=0.0
               fa(3*(ik-1)+1)=0.0
               fa(3*(ik-1)+2)=0.0
               fa(3*(ik-1)+3)=0.0

cd         write(lunout,*)  '  FPI endatoms:  nFPI = ',nFPI
cd           write(lunout,*) '  After rigid atom operation: time=',time
cd           write(lunout,233) 1,(va(3*k),k=1,natms(1))
cd           write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))

            endif
c     close `if(itagFPI(ip) .ge. 100....) then'            

c     Chain FPI with constrained centroid:
                                                 
            if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then
c     for tags 300 to 399 zero the net velocity and force on the
                   
c     chain of images as a whole:
                                                  
               vxtot=0.0
               vytot=0.0
               vztot=0.0
cc           fxtot=0.0
cc           fytot=0.0
cc           fztot=0.0

               fact=1./nimFPI(ip)
               do iim=1,nimFPI(ip)
                  ind=3*(iFPI(ip)+(iim-1)-1)
c     add up velocity of all images in this chain:
                  vxtot=vxtot+va(ind+1)
                  vytot=vytot+va(ind+2)
                  vztot=vztot+va(ind+3)
c     add up force on all images in this chain:
cc             fxtot=fxtot+fa(ind+1)
cc             fytot=fytot+fa(ind+2)
cc             fztot=fztot+fa(ind+3)
               enddo
c     now subtract the total from the vecor for each image:
               do iim=1,nimFPI(ip)
                  ind=3*(iFPI(ip)+(iim-1)-1)
c     subtract the total velocity from vel. of all images in this chain:
c$$$                  v2perm=v2perm+va(3*ind+1)**2+va(3*ind+2)**2+
c$$$     +                 va(3*ind+3)**2
c     'v2perm' keeps track of the total kinetic energy taken out

                  va(ind+1)=va(ind+1)-vxtot*fact
                  va(ind+2)=va(ind+2)-vytot*fact
                  va(ind+3)=va(ind+3)-vztot*fact
                   
c$$$                  v2perm=v2perm-(va(ind+1)**2+va(ind+2)**2+va(ind+3)**2)

c     subtract the total force from force on all images in this chain:
cc              fa(ind+1)=fa(ind+1)-fxtot
cc              fa(ind+2)=fa(ind+2)-fytot
cc              fa(ind+3)=fa(ind+3)-fztot
               enddo

            endif
c     close `if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then'

 137     continue

      endif                     ! close `if(nFPI .gt. 0) then'

c-----------------------------------------------------------------------
c  Constrained atoms:  tags 4, 5, 6, 7, 8 and 9
c         and images in FPI with tags 540-599.

      chainLength = nimFPI(1)
      if (nimFPI(1) .eq. 0) chainlength = 1
      totM = 2.d0 * aMass(1) + aMass(2)
      do 140 ip=1,naconstr
         ik=iconstr(ip)
         if((itag(ik) .eq. 4) .or.
     +        (itag(ik).ge.540.and.itag(ik).le.549))then
c     zero the y and z components
c$$$            v2perm=v2perm+va(3*(ik-1)+2)**2+va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+2)=0.0
            va(3*(ik-1)+3)=0.0
            fa(3*(ik-1)+2)=0.0
            fa(3*(ik-1)+3)=0.0
         endif
         if((itag(ik) .eq. 5) .or.
     +        (itag(ik).ge.550.and.itag(ik).le.559))then
c     zero the x and z components
c$$$            v2perm=v2perm+va(3*(ik-1)+1)**2+va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+1)=0.0
            va(3*(ik-1)+3)=0.0
            fa(3*(ik-1)+1)=0.0
            fa(3*(ik-1)+3)=0.0
         endif
         if((itag(ik) .eq. 6) .or.
     +        (itag(ik).ge.560.and.itag(ik).le.569))then
c     zero the x and y components
c$$$            v2perm=v2perm+va(3*(ik-1)+1)**2+va(3*(ik-1)+2)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+1)=0.0
            va(3*(ik-1)+2)=0.0
            fa(3*(ik-1)+1)=0.0
            fa(3*(ik-1)+2)=0.0
         endif
         if((itag(ik) .eq. 7) .or.
     +        (itag(ik).ge.570.and.itag(ik).le.579))then
c     zero the xcomponent
c$$$            v2perm=v2perm+va(3*(ik-1)+1)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+1)=0.0
            fa(3*(ik-1)+1)=0.0
         endif
         if((itag(ik) .eq. 8) .or.
     +        (itag(ik).ge.580.and.itag(ik).le.589))then
c     zero the xcomponent
c$$$            v2perm=v2perm+va(3*(ik-1)+2)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+2)=0.0
            fa(3*(ik-1)+2)=0.0
         endif
         if((itag(ik) .eq. 9) .or.
     +        (itag(ik).ge.590.and.itag(ik).le.599))then
c     zero the xcomponent
c$$$            v2perm=v2perm+va(3*(ik-1)+3)**2
c     'v2perm' keeps track of the total kinetic energy taken
            va(3*(ik-1)+3)=0.0
            fa(3*(ik-1)+3)=0.0
         endif

cERB Molecules with fixed center of mass.
         if(itag(ik) .eq. 11) then
            ik = ik - nAtms(1)
            moleNumber = int((ik-1)/chainLength) + 1
            image      = ik - chainLength * (moleNumber - 1)
            ik = moleNumber

            indexO = 3 * (nAtms(1) + (ik-1) * chainLength + 
     $           (image-1))
            indexH1 = 3 * (2 * (ik-1) * chainLength + (image-1))
            indexH2 = indexH1 + 3 * chainLength
            do i = 1, 3
               vCM(i) = (amass(1) * (va(indexH1+i) + va(indexH2+i)) + 
     $              amass(2) * va(indexO+i)) / totM
               fTotal(i) = fa(indexH1+i) + fa(indexH2+i) + fa(indexO+i)
            end do
            do i = 1, 3
               va(indexH1+i) = va(indexH1+i) - vCM(i)
               va(indexH2+i) = va(indexH2+i) - vCM(i)
               va(indexO+i) = va(indexO+i) - vCM(i)

               fa(indexH1+i) = fa(indexH1+i) - fTotal(i) * amass(1)/totM
               fa(indexH2+i) = fa(indexH2+i) - fTotal(i) * amass(1)/totM
               fa(indexO+i) = fa(indexO+i) - fTotal(i) * amass(2)/totM
            end do
         end if

cERB Molecule with its CM constrained to a plane of normal n.
         if((itag(ik) .ge. 12) .and. (itag(ik) .le. 16)) then
            if (itag(ik) .eq. 12) then
c     CM constrained to y-z plane (x fixed)
               n(1) = 1.d0
               n(2) = 0.d0
               n(3) = 0.d0
            else if(itag(ik) .eq. 13) then
c     CM constrained to x-z plane (y fixed)
               n(1) = 0.d0
               n(2) = 1.d0
               n(3) = 0.d0
            else if(itag(ik) .eq. 14) then
c     CM constrained to x-y plane (z fixed)
               n(1) = 0.d0
               n(2) = 0.d0
               n(3) = 1.d0
            else if(itag(ik) .eq. 15) then
c     CM constrained to the vertical plane at 30 degrees with x axis
               n(1) = -.5d0
               n(2) = 0.866025403784438597d0
               n(3) = 0.d0
            else if(itag(ik) .eq. 16) then
c     CM constrained to the vertical plane at 150 degrees with x axis
               n(1) = .5d0
               n(2) = 0.866025403784438597d0
               n(3) = 0.d0
            end if

            ik = ik - nAtms(1)
            moleNumber = int((ik-1)/chainLength) + 1
            image      = ik - chainLength * (moleNumber - 1)
            ik = moleNumber

            indexO = 3 * (nAtms(1) + (ik-1) * chainLength + 
     $           (image-1))
            indexH1 = 3 * (2 * (ik-1) * chainLength + (image-1))
            indexH2 = indexH1 + 3 * chainLength
            do i = 1, 3
               vCM(i) = (amass(1) * (va(indexH1+i) + va(indexH2+i)) + 
     $              amass(2) * va(indexO+i)) / totM
               fTotal(i) = fa(indexH1+i) + fa(indexH2+i) + fa(indexO+i)
            end do

            vDotN = vCM(1)*n(1) + vCM(2)*n(2) + vCM(3)*n(3)
            fDotN = fTotal(1)*n(1) + fTotal(2)*n(2) + fTotal(3)*n(3)
c     Subtract the perpendicular component of the force and velocity of
c     the CM.
            do i = 1, 3
               vCM(i) = n(i) * vDotN
               fTotal(i) = n(i) * fDotN
            end do
            do i = 1, 3
               va(indexH1+i) = va(indexH1+i) - vCM(i)
               va(indexH2+i) = va(indexH2+i) - vCM(i)
               va(indexO+i) = va(indexO+i) - vCM(i)

               fa(indexH1+i) = fa(indexH1+i) - fTotal(i) * amass(1)/totM
               fa(indexH2+i) = fa(indexH2+i) - fTotal(i) * amass(1)/totM
               fa(indexO+i) = fa(indexO+i) - fTotal(i) * amass(2)/totM
            end do
         end if

cERB Molecule allowed to move just along the z axis.
         if((itag(ik) .ge. 17) .and. (itag(ik) .le. 19)) then
            if (itag(ik) .eq. 17) then
c     CM constrained to x axis.
               n(1) = 0.d0
               n(2) = 1.d0
               n(3) = 1.d0
            else if(itag(ik) .eq. 18) then
c     CM constrained to y axis
               n(1) = 1.d0
               n(2) = 0.d0
               n(3) = 1.d0
            else if(itag(ik) .eq. 19) then
c     CM constrained to z axis
               n(1) = 1.d0
               n(2) = 1.d0
               n(3) = 0.d0
            end if
            ik = ik - nAtms(1)
            moleNumber = int((ik-1)/chainLength) + 1
            image      = ik - chainLength * (moleNumber - 1)
            ik = moleNumber

            indexO = 3 * (nAtms(1) + (ik-1) * chainLength + 
     $           (image-1))
            indexH1 = 3 * (2 * (ik-1) * chainLength + (image-1))
            indexH2 = indexH1 + 3 * chainLength
            do i = 1, 3
               vCM(i) = (amass(1) * (va(indexH1+i) + va(indexH2+i)) + 
     $              amass(2) * va(indexO+i)) / totM * n(i)
               fTotal(i) = (fa(indexH1+i) + fa(indexH2+i) + 
     $              fa(indexO+i)) * n(i)
            end do
            do i = 1, 3
               va(indexH1+i) = va(indexH1+i) - vCM(i)
               va(indexH2+i) = va(indexH2+i) - vCM(i)
               va(indexO+i) = va(indexO+i) - vCM(i)

               fa(indexH1+i) = fa(indexH1+i) - fTotal(i) * amass(1)/totM
               fa(indexH2+i) = fa(indexH2+i) - fTotal(i) * amass(1)/totM
               fa(indexO+i) = fa(indexO+i) - fTotal(i) * amass(2)/totM
            end do
         end if
 140  continue

c-----------------------------------------------------------------------
      if (.false.) then
      oi = 1
      totM = 2.d0 * aMass(1) + aMass(2)
      indexH1 = 3 * (2 * oi - 2)
      indexH2 = 3 * (2 * oi - 1)
      indexO = 3 * (16 + oi - 1)
      do i = 1, 3
         rCMi(i) = (amass(1) * (ra(indexH1+i) + ra(indexH2+i)) + 
     $        amass(2) * ra(indexO+i)) / totM
      end do
      do k = 2, 4
         indexH1 = 3 * (2 * k - 2)
         indexH2 = 3 * (2 * k - 1)
         indexO = 3 * (16 + k - 1)
      
         do i = 1, 3
            rCMk(i) = (amass(1) * (ra(indexH1+i) + ra(indexH2+i)) + 
     $           amass(2) * ra(indexO+i)) / totM
            vCM(i) = (amass(1) * (va(indexH1+i) + va(indexH2+i)) + 
     $           amass(2) * va(indexO+i)) / totM
            fTotal(i) = fa(indexH1+i) + fa(indexH2+i) + fa(indexO+i)
         end do
         
         do i = 1, 3
            dr(i) = rCMk(i) - rCMi(i)
         end do

         ddd = sqrt(dot(dr, dr))
         do i = 1, 3
            dr(i) = dr(i) / ddd
         end do
         fTddr = dot(ft, dr)
         vDotN = dot(vCM, dr)
         do i = 1, 3
            fTotal(i) = fTotal(i) - fTddr * dr(i)
            vCM(i) = vCM(i) - vDotN * dr(i)
         end do

c     Subtract the perpendicular component of the force and velocity of
c     the CM.
         do i = 1, 3
            va(indexH1+i) = va(indexH1+i) - vCM(i) * totM / amass(1)/ 3.
            va(indexH2+i) = va(indexH2+i) - vCM(i) * totM / amass(1)/3.
            va(indexO+i) = va(indexO+i) - vCM(i) * totM / amass(2)/3.
            
            fa(indexH1+i) = fa(indexH1+i) - fTotal(i) * amass(1)/totM
            fa(indexH2+i) = fa(indexH2+i) - fTotal(i) * amass(1)/totM
            fa(indexO+i) = fa(indexO+i) - fTotal(i) * amass(2)/totM
         end do
      end do

c--------------------------------------------------

      oi = 5
      indexH1 = 3 * (2 * oi - 2)
      indexH2 = 3 * (2 * oi - 1)
      indexO = 3 * (16 + oi - 1)
      do i = 1, 3
         rCMi(i) = (amass(1) * (ra(indexH1+i) + ra(indexH2+i)) + 
     $        amass(2) * ra(indexO+i)) / totM
      end do
      do k = 6, 8
         indexH1 = 3 * (2 * k - 2)
         indexH2 = 3 * (2 * k - 1)
         indexO = 3 * (16 + k - 1)
      
         do i = 1, 3
            rCMk(i) = (amass(1) * (ra(indexH1+i) + ra(indexH2+i)) + 
     $           amass(2) * ra(indexO+i)) / totM
            vCM(i) = (amass(1) * (va(indexH1+i) + va(indexH2+i)) + 
     $           amass(2) * va(indexO+i)) / totM
            fTotal(i) = fa(indexH1+i) + fa(indexH2+i) + fa(indexO+i)
         end do
         
         do i = 1, 3
            dr(i) = rCMk(i) - rCMi(i)
         end do

         ddd = sqrt(dot(dr, dr))
         do i = 1, 3
            dr(i) = dr(i) / ddd
         end do
         fTddr = dot(ft, dr)
         vDotN = dot(vCM, dr)
         do i = 1, 3
            fTotal(i) = fTotal(i) - fTddr * dr(i)
            vCM(i) = vCM(i) - vDotN * dr(i)
         end do

c     Subtract the perpendicular component of the force and velocity of
c     the CM.
         do i = 1, 3
            va(indexH1+i) = va(indexH1+i) - vCM(i) * totM / amass(1)/ 3.
            va(indexH2+i) = va(indexH2+i) - vCM(i) * totM / amass(1)/3.
            va(indexO+i) = va(indexO+i) - vCM(i) * totM / amass(2)/3.
            
            fa(indexH1+i) = fa(indexH1+i) - fTotal(i) * amass(1)/totM
            fa(indexH2+i) = fa(indexH2+i) - fTotal(i) * amass(1)/totM
            fa(indexO+i) = fa(indexO+i) - fTotal(i) * amass(2)/totM
         end do
      end do

      end if
      
      return
      end
