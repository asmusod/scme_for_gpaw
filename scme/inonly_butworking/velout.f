c   EAM e93:
c   EAM version c92:  reserve tags 100-199 for FPI with fixed end atoms.                  
c                        -     -   200-299  -   -  that form closed loop.                 
c                        -     -   300-399  -   - loops with constrained centroid         
c   EAM version b91: more atoms allowed and add FPI with fixed ends (tag 10-19)           
c                                                                                         
C THIS SUBROUTINE WRITES OUT THE CONFIGURATION FILE:                                      
c   indout:    This parameter determines what gets written out.                           
c    indout < 0  -  the velocities and potential are not written out.           
c    |indout| = 2 then raNoPBC gets written, not RA                  
c                 (raNoPBC are the coordinates without periodic boundary conditions)            
c                 and the potential energy of each atom gets written out in the                
c                 fourth column after the velocities. (this is used in CM.CON files)           
c    |indout| = 3  then only atoms with z>zminwr get written out. 
C                                                                                         
c    cooscale:    coordinates get multiplied by cooscale before writing                   

      SUBROUTINE velOUT(lun, cooscale, iVelAtms, nMolSave)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comoutcntrl.cmn'

      COMMON /ISEEDS/ ISEED1

      dimension ndumpatms(MAXCOMP), iVelAtms(*)

c    -------------------------------------------------------------------------

c    Debug write statements:                                                              

      WRITE(lun,100) ISEED1
 100  FORMAT(i20,20x,'                       RANDOM NUMBER SEED')
C
      write(lun,101) time
 101  format(f14.4,20x,'                               Time')
      write(lun,102) ax,ay,az
 102  format(g23.15,g23.15,g23.15,'ax,ay,az')
      write(lun,104) alpha,beta,gamma
 104  format(g16.8,g16.8,g16.8,10x,'     alpha,beta,gamma')
      write(lun,105) pmass,volv
 105  format(g24.15,g24.15,10x,'     Pmass, Volv')
      write(lun,1003) naperm,nharm,nFPI
 1003 format(i5,2x,i5,2x,i5, '       naperm,nharm,nFPI:')
      WRITE(lun,110) NATYPE
 110  FORMAT(i10,45x,'     # OF TYPES OF ATOMS')
                                                                                       
c   The usual case:
      WRITE(lun,120) (NATMS(I) * nMolSave/natms(2),I=1,NATYPE)
 120  FORMAT(10i5,45x,'     # OF ATOMS OF EACH TYPE')

      WRITE(lun,130) (AMASS(I),I=1,NATYPE)
 130  FORMAT(10f8.3,' MASS(i)')
c                                                                                         
      isumat=0

      i = 1
      WRITE(lun,156)(CPHEAD(I,K),K=1,78)
 156  FORMAT(1X,80A1)
      WRITE(lun,157) I
 157  FORMAT('   COORDINATES OF COMPONENT #',I2)
         
      do j = 1, nMolSave
         l = iVelAtms(j)
         k = 2*(l-1) + 1
         x = VA(3*K-2) * cooscale
         y = VA(3*K-1) * cooscale
         z = VA(3*K) * cooscale
         WRITE(lun,140) x,y,z,itag(K),K

         k = k + 1
         x = VA(3*K-2) * cooscale
         y = VA(3*K-1) * cooscale
         z = VA(3*K) * cooscale
         WRITE(lun,140) x,y,z,itag(K),K
      end do

      i = 2
      WRITE(lun,156)(CPHEAD(I,K),K=1,78)
      WRITE(lun,157) I
         
      do j = 1, nMolSave
         l = iVelAtms(j)
         k = l + nAtms(1)
         x = VA(3*K-2) * cooscale
         y = VA(3*K-1) * cooscale
         z = VA(3*K) * cooscale
         WRITE(lun,140) x,y,z,itag(K),K
      end do


 140     FORMAT(3g23.15,i5,i6)

      return

      END
