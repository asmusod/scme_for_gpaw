C     THIS SUBROUTINE COMPUTES THE KINETIC ENERGY
c     The results appear on the commonblock 'ENERGY'.
C     On exit the variable TTOT contains the sum of   MASS*(VELOCITY)**2
c     NB!  the velocities are scaled , so the kinetic energy is multiplied
c     by ax**2. 
c     Also, unlike older versions of KINET, this version multiplies by 0.5
c     so the output is really the kinetic energy.
c     UT(1) is potential energy of atoms of type 1 interacting with themselves
c     UT(3)  -    -        -          -       -  2      -        -      - 
c     UT(2)  -    -        -          -       -  1 and type 2 interacting

      SUBROUTINE KINET

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comtgr.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comoutcntrl.cmn'
      include '../commonblks/comluns.cmn'


cPres:   scale the mesh sizes:

      tx=tgrlx/ax
      ty=tgrly/ax
      tz=tgrlz/ax

cc      write(6,*) '    KINET:     time = ',time,'  VA = ',(va(k),k=1,6)

      nuFPI = nimFPI(1)
      if (nuFPI .eq. 0) nuFPI = 1
      do i = 1, nuFPI
         T(i) = 0.d0
      end do

      INDX=0
      DO 100 ITYPE=1,NATYPE
         fact = 0.5d0 * AMASS(ITYPE) * ax*ax

         DO Ia = 1, NATMS(ITYPE)/nuFPI
            do jFPI = 1, nuFPI

               ind = 3 * nuFPI * (ia-1) + 3 * (jFPI-1) + indx
               vx = VA(1+ind)
               vy = VA(2+ind)
               vz = VA(3+ind)
               T(jFPI) = T(jFPI) + fact * (vx*vx + vy*vy + vz*vz)

            end do              ! close jFPI
         end do                 ! close loop over atoms of component 'ITYPE'.

         INDX=INDX+3*NATMS(ITYPE)
 100  CONTINUE
      
      TTOT = 0.d0
      do jFPI = 1, nuFPI
         TTOT = TTOT + T(jFPI)
      end do
c     

      RETURN

 995  continue
      write(lunout,200) ix,iy,iz,tgrlx,tgrly,tgrlz,(indx/3+ia),
     +             ra(1+indx+3*(ia-1)),ra(2+indx+3*(ia-1)),
     +             ra(3+indx+3*(ia-1))
200   format(/' ERROR 995 kinet: mesh size too small, ix,iy,iz=',3i5/
     +        ' tgrlx,tgrly,tgrlz = ',3g12.4,'  atom i=',i5,' has '/
     +        '   x,y,z = ',3g15.5)
      return
      END
