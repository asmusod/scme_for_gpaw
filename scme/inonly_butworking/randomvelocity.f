c                                                                                         
C  THIS ROUTINE GENERATES RANDOM VELOCITIES FOR AN ATOM OF A GIVEN MASS, 'AMASS'          
C  CHARACTERISTIC OF A GIVEN BATH TEMPERATURE, 'TEMP'.                                    
c   The velocity components are written into the array 'V'.                               
c   NB:  the velocities are not scaled.                                                   

      SUBROUTINE RANVEL(Vx,Vy,Vz,AMASS,TEMP)
c                                                                                         
      implicit real*8 (a-h,o-z)
c                                                                                         
c      include '../commonblks/parameters.cmn'
c      include '../commonblks/comgeom.cmn'
c      include '../commonblks/comluns.cmn'

      COMMON/ISEEDS/ISEED1
c      data convFactor / 1.018050697d0 /
      data convFactor / 1.d0 /


      if(temp .lt. 0.1e-10) then
           Vx=0.0
           Vy=0.0
           Vz=0.0
           return
      endif

c      ftemp = 2./3.
c----------------------- ERB
      ftemp = 1.d0
c-----------------------
      if(az .le. 0.99e-06) ftemp= 1.0

c                                                                                         
cq      write(lunout,221)                                                                 
cq221   format(/'      ranvel:  '/                                                        
cq     +        '  iseedx:  ranx:    iseedy:  rany:    iseedz:  ranz:')                   
c                                                                                         
        FACT=SQRT(2.0*(TEMP*ftemp)/AMASS)

        ranx=RANb(ISEED1,ISEED2)
c        Vx=FACT*ERFINV(2.*ranx-1.) * convFactor
        Vx=FACT*ERFINV(2.*ranx-1.)
        rany=RANb(ISEED1,ISEED2)
c        Vy=FACT*ERFINV(2.*rany-1.) * convFactor
        Vy=FACT*ERFINV(2.*rany-1.)
        ranz=RANb(ISEED1,ISEED2)
c        Vz=FACT*ERFINV(2.*ranz-1.) * convFactor
        Vz=FACT*ERFINV(2.*ranz-1.)
c                                                                                         
cd        write(lunout,220) vx,vy,vz                                                      
cd220     format('   RANVEL:    vx,vy,vz= ',g12.3,'   ',g12.3,'  ',g12.3)                 
c                                                                                         
       RETURN
       END

c-----------------------------------------------------------------------

      SUBROUTINE RANw(w, J, TEMP)

c     Get three components for the angular velocity of the rigid H2O
c     molecule. Each component is drawn at random from a gaussian
c     distribution exp(-1/2 J w*w/kT) where J is the moment of inertia
c     respect to each principal axis.

      implicit real*8 (a-h,o-z)

c      include '../commonblks/parameters.cmn'
c      include '../commonblks/comgeom.cmn'
c      include '../commonblks/comluns.cmn'

      COMMON/ISEEDS/ISEED1
      data convFactor / 1.018050697d0 /

      real*8 J(3), w(3), ranx
c$$$      print '(A)', 'Inertia moments'
c$$$      print '(3f15.10)', j(1), j(2), j(3)
c$$$      pause ' at ranW'

      if(temp .lt. 0.1e-10) then
         do i = 1, 3
            w(i) = 0.d0
         end do
      else
c      ftemp = 2./3.
c----------------------- ERB
         ftemp = 1.d0
c-----------------------
         do i = 1, 3
            FACT = SQRT(2.d0 * (TEMP * ftemp) / J(i))
            ranx = RANb(ISEED1,ISEED2)
            w(i) = FACT * ERFINV(2. * ranx - 1.d0)
         end do
      end if

      RETURN
      END
