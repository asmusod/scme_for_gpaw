c   Centroid constraints:                                                                          
c      this subroutine zeros the force on the centroid of FPI chains with
c      tags between 300 and 399.                              

       subroutine centroidconstr

       implicit real*8 (a-h,o-z)

       include '../commonblks/parameters.cmn'
       include '../commonblks/comgeom.cmn'
       include '../commonblks/comconf.cmn'
       include '../commonblks/comtime.cmn'
       include '../commonblks/comluns.cmn'

cd      write(6,*) '  CENTROID:  time = ',time,'  naperm = ',naperm,                         
cd     +       '   iperm = ',(iperm(k),k=1,naperm)                                        

c     -----------------------------------------------------------------                   
c     -----------------------------------------------------------------                   

c  FPI chains:                                                                            
      if(nFPIcnt .gt. 0) then
         do 137 ip=1,nFPI
c          loop over FPI chains:                                                          

           if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then
c            Here is a chain FPI with constrained centroid:                                                 
c            for tags 300 to 399 zero the net force on the                   
c            chain of images as a whole:                                                  
           fxtot=0.0
           fytot=0.0
           fztot=0.0
           fact=1./nimFPI(ip)
           do iim=1,nimFPI(ip)
             ind=3*(iFPI(ip)+(iim-1)-1)
c              add up force on all images in this chain:                                     
             fxtot=fxtot+fa(ind+1)
             fytot=fytot+fa(ind+2)
             fztot=fztot+fa(ind+3)
           enddo
c           now subtract the total from the vector for each image:                         
           do iim=1,nimFPI(ip)
             ind=3*(iFPI(ip)+(iim-1)-1)
c              subtract the total velocity from vel. of all images in this chain:            
cv             v2perm=v2perm+va(3*ind+1)**2+va(3*ind+2)**2+
cv     +                    va(3*ind+3)**2
c                'v2perm' keeps track of the total kinetic energy taken out               

cv           va(ind+1)=va(ind+1)-vxtot
cv           va(ind+2)=va(ind+2)-vytot
cv           va(ind+3)=va(ind+3)-vztot

cv             v2perm=v2perm-(va(ind+1)**2+va(ind+2)**2+va(ind+3)**2)

c           subtract the total force from force on all images in this chain:              
              fa(ind+1)=fa(ind+1)-fxtot*fact
              fa(ind+2)=fa(ind+2)-fytot*fact
              fa(ind+3)=fa(ind+3)-fztot*fact
           enddo

c      Store the total force on the FPI chain in the array cntrFA:
           cntrFA(3*(ip-1)+1)=fxtot
           cntrFA(3*(ip-1)+2)=fytot
           cntrFA(3*(ip-1)+3)=fztot

          endif
c             close `if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then'            


137       continue
c                     close loop over FPI chains
          endif
c             close `if(nFPI .gt. 0) then'                                                

c     -----------------------------------------------------------------                   

        return
      end

