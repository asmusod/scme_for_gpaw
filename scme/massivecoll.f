c   This subroutine performs massive stochastic collisions or zeroes all                  
c     velocities, depending on temperature TEXTL.                                         

      Subroutine masscol(textl,pkinet,ekindiff)

      implicit real*8(a-h,o-z)
                                                                                        
      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comluns.cmn'

      COMMON / REFENC / REFPRE,REFVOL,REFPOT,REFKIN,REFENT,REFENG


      ftemp = 2./3.
      if(az .le. 0.99e-06) ftemp=1.0
c                     ftemp=1 is appropriate for 2-D simulations.                         

      axinv=1./ax
c      write(*,*)'tempextl = ',textl                                                      

cq         write(lunout,300) ax,ay,az,textl,pkinet                                        
cq300      format(/'   masscol:  ax,ay,az,textl,pkinet=',5g12.4)                          
c                                                                                         
      CALL KINET
c              calculate the total kinetic energy before massive collision for            
c              bookkeeping.                                                               
      ekinbefore=ttot
c                                                                                         
          ISTR=1
          if(textl .gt. 0.1e-04) then
            DO 210 ITYPE=1,NATYPE
              amasst=AMASS(ITYPE)
              do iat=1,natms(itype)
                 ind=ISTR+3*(iat-1)
                 CALL RANVEL(VA(ind),VA(ind+1),VA(ind+2),AMASSt,TEXTL)
              enddo
              ISTR=ISTR+3*NATMS(ITYPE)
210           continue
cPres:      scale all the velocities:                                                     
              do icoo=1,3*natoms
                  va(icoo)=va(icoo)*axinv
              enddo
          else
c           Zero the velocity of all atoms (steepest descent):                            
            do 230 i=1,natoms
               va(3*(i-1)+1)=0.0
               va(3*(i-1)+2)=0.0
               va(3*(i-1)+3)=0.0
230            continue
          endif
c                                                                                         
c      write info:                                                                        
       if(textl .gt. 0.001) write(lunout,215)  ftemp*textl,time
215    format('#   Massive collision with',
     +        '  T bath = ',g12.4,' at time =',f10.3)
c                                                                                         
cq       sumva2=0.                                                                        
cq       do 217 i=1,natoms                                                                
cq         va2=va(3*i-2)**2+va(3*i-1)**2+va(3*i)**2                                       
cq         sumva2=sumva2+va2                                                              
cq         write(lunout,216) i,(va(3*i-2),va(3*i-1),va(3*i)),va2                          
cq216      format('   ',i5,4g12.4)                                                        
cq217      continue                                                                       
cq       write(lunout,218) sumva2                                                         
cq218    format(/'        sumVa2 = ',g12.4/)                                              
c                                                                                         
c           --              --                --             --                           
cd      write(lunout,233) 'A Vz',(va(3*(k)),k=1,natms(1))                                 
cd      write(lunout,233) 'B Vz',(va(3*(natms(1)+k)),k=1,natms(2))                        
cd233   format('   Masscol:  Vz of ',a4,' atoms:'/                                        
cd     +        5g15.6/5g15.6/5g15.6/5g15.6)                                              
cd      write(lunout,234) naperm,(iperm(i),i=1,naperm)                                    
cd234   format('   naperm = ',i4,'  iperm = ',14i5/14i5/14i5)                             
c                                                                                         
c     Rigid:                                                                              
          call rigidatoms(dummy)
c               zero the velocity of rigid (permafrost) atoms:                            

c     MD2D operation:                                                                     
          if(az .le. 0.99e-06) call md2dop(dummy)
c                zero the velocity in z direction if this is a 2D simulation.             
c                                                                                         
cd      write(lunout,233) 'A Vz',(va(3*(k)),k=1,natms(1))                                 
cd      write(lunout,233) 'B Vz',(va(3*(natms(1)+k)),k=1,natms(2))                        
c            --             --                --              --                          
c                                                                                         
      CALL KINET
c                                                                                         
      ekinafter=ttot
      ekindiff=ekinafter-ekinbefore
c                                                                                         
      PKNEW=TTOT
      pkinchng=(PKINET-PKNEW)
cq      RANKIN=RANKIN+pkinchng                                                            
c          this is the change in kinetic energy due to the massive collision.             
cq      PKINET=PKNEW                                                                      
c                                                                                         
       pkinet=ekinafter
c                                                                                         
cq       write(lunout,219)  ekindiff,pkinet                                               
cq219    format('   ekindiff = ',g20.10,'  pkinet = ',g20.10)                             
c                                                                                         
ccq       if(dabs(ekindiff+pkinchng) .gt. 0.1e-5) go to 999                               
c                                                                                         
       return
c                                                                                         
ccq999    continue                                                                        
ccq       write(lunout,9999) ekindiff,pkinchng                                            
ccq9999   format(/'   ERROR masscoll:  ekindiff = ',g20.10,'  pkinchng = ',               
ccq     +             g20.10)                                                             
ccq      return                                                                           
c                                                                                         
      end


c    hit end in read 100, ic =           112                                              

c      total number of lines     =        117                                             
c      number of new lines       =          0                                             
c      number of type *, found   =          0                                             
 
 

c   Output from fortfixIBM (fixes line lengths):
c      hit end in read 100,   ic =        121
c      total number of lines     =        123
c      number of new lines       =          0
c      number of tabs found      =          0
c      number of type *, found   =          0
c-----------------------------------------------------------------------

c   This subroutine, depending on temperature TEXTL, performs massive
c   stochastic collisions or zeroes all velocities, of the rigid water
c    molecules  

      Subroutine masscolH2O(textl, pkinet, ekindiff, InertiaMoment, 
     $     mMolecule)
      
      implicit real*8(a-h,o-z)
                                                                                        
      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comluns.cmn'

      COMMON / REFENC / REFPRE,REFVOL,REFPOT,REFKIN,REFENT,REFENG

      real*8 r1(3), r2(3), r3(3), v1(3), v2(3), v3(3)
      real*8 InertiaMoment(3), angularVel(3), vCM(3), mMolecule
      integer chainLength

c$$$      print '(A)', 'Inertia moments'
c$$$      print '(3f15.10)', InertiaMoment(1), InertiaMoment(2), 
c$$$     $     InertiaMoment(3)
c$$$
c$$$      print '(A, f10.5)', 'molecule mass = ', mMolecule
c$$$      print '(A, f10.5)', 'temperature = ', textl

      axinv=1./ax
      chainLength = nimFPI(1)
      if (nimFPI(1) .eq. 0) chainlength = 1

c     calculate the total kinetic energy before massive collision for
c     bookkeeping.                     

      CALL KINET
      ekinbefore=ttot

      if(textl .gt. 0.1e-04) then


         iseed1 = 1
         iseed2 = 987
         
         nOxygens = nAtms(2)/chainLength
         nHydrogens = natms(1)/chainLength

         do iatom = 1, nOxygens

            do image = 1, chainLength

               indexO = 3 * (nAtms(1) + (iatom-1) * chainLength + 
     $              (image-1))
               indexH1 = 3 * (2 * (iatom-1) * chainLength + (image-1))
               indexH2 = indexH1 + 3 * chainLength
               do ii = 1, 3
                  r1(ii) = RA(indexH1 + ii)
                  r2(ii) = RA(indexH2 + ii)
                  r3(ii) = RA(indexO + ii)
               end do

c            print *, ' Colliding massively with the bath at ', textl
c            pause 'at masscloH2O'

c     Draw a random velocity for the center of mass of the molecule
               CALL RANVEL(vCM(1), vCM(2), vCM(3), mMolecule, TEXTL)

c     Scale the velocity of the center of mass dividing by ax
               do ii = 1, 3
                  vCM(ii) = vCM(ii) * axinv
c$$$               vCM(ii) = 0.d0
               end do

c$$$            print '(3(f15.10,2x))', (vCM(j), j=1,3)
c$$$            print *

c     Get a random velocity for the angular velocity of the molecule, in
c     reference frame of its principal axes.
               CALL RANw(angularVel, InertiaMoment, TEXTL)

c$$$            print '(3(f15.10,2x))', (angularVel(j), j=1,3)
c$$$            pause 'at masscolH2O'

c$$$            angularVel(1) = 0.d0
c$$$            angularVel(2) = 0.d0
c$$$            angularVel(3) = 0.d0

c     Compute the tangential velocity of each atom (r x w) in the
c     reference frame fixed to the crystal.                      
               call tangentialVel(r1, r2, r3, v1, v2, v3, angularVel, 
     $              mMolecule)

c$$$            do ii= 1, 3
c$$$               v1(ii) = 0.d0
c$$$               v2(ii) = 0.d0
c$$$               v3(ii) = 0.d0
c$$$            end do

c     Store the new velocities in VA, V = Vcm + Vtangential
               indexO = 3 * (nAtms(1) + (iatom-1) * chainLength + 
     $              (image-1))
               indexH1 = 3 * (2 * (iatom-1) * chainLength + (image-1))
               indexH2 = indexH1 + 3 * chainLength
               do ii = 1, 3
                  VA(indexH1 + ii) = v1(ii) + vCM(ii)
                  VA(indexH2 + ii) = v2(ii) + vCM(ii)
                  VA(indexO + ii)  = v3(ii) + vCM(ii)
               end do
            end do
         end do
      else 
c     Zero the velocity of all atoms, in case we are doing steepest descent:
         do i = 1, natoms
            va(3*(i-1)+1)=0.0
            va(3*(i-1)+2)=0.0
            va(3*(i-1)+3)=0.0
         end do
      endif

      if(textl .gt. 0.001) write(lunout,215)  textl,time
 215  format('#   Massive collision with',
     +     '  T bath = ',g12.4,' at time =',f10.3)

      call rigidatoms(dummy)
      CALL KINET

      ekinafter=ttot
      ekindiff=ekinafter-ekinbefore

      PKNEW=TTOT
      pkinchng=(PKINET-PKNEW)

      pkinet=ekinafter

c$$$      write(*,219)  ekindiff,pkinet
c$$$ 219  format('   ekindiff = ',g20.10,'  pkinet = ',g20.10)

      return
      end







