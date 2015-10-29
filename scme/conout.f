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

      SUBROUTINE CONOUT(lun,indout,cooscale)

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comoutcntrl.cmn'

      COMMON /ISEEDS/ ISEED1

      dimension ndumpatms(MAXCOMP)

c    -------------------------------------------------------------------------

c    Debug write statements:                                                              
cd      write(6,*)                                                                        
cd      write(6,*) ' ** CONOUT:  Natoms =',Natoms,' Cooscale =',cooscale                  
cd      do i=1,Natoms                                                                     
cd        write(6,*) (RA(3*(i-1)+k),k=1,3)                                                
cd      enddo                                                                             
cd      write(6,*)                                                                        

      if(iseed1 .eq. 0) then
         WRITE(lun,121)
121      FORMAT(1x,'  10000 ',30x,'        default  RANDOM NUMBER SEED')
      else
         WRITE(lun,100) ISEED1
100      FORMAT(i20,20x,'                       RANDOM NUMBER SEED')
      endif
C                                                                                         
      write(lun,101) time
101   format(f14.4,20x,'                               Time')
      write(lun,102) ax,ay,az
102   format(g23.15,g23.15,g23.15,'ax,ay,az')
      write(lun,104) alpha,beta,gamma
104   format(g16.8,g16.8,g16.8,10x,'     alpha,beta,gamma')
      write(lun,105) pmass,volv
105   format(g24.15,g24.15,10x,'     Pmass, Volv')
cERB      write(lun,103) naperm,nharm,nFPI,(iFPI(if),if=1,nFPI)
cERB 103   format(i5,2x,i5,2x,i5,300i4,
cERB     +                   '       naperm,nharm,nFPI:')
      write(lun,1003) naperm,nharm,nFPI
 1003 format(i5,2x,i5,2x,i5, '       naperm,nharm,nFPI:')
      WRITE(lun,110) NATYPE
110   FORMAT(i10,45x,'     # OF TYPES OF ATOMS')
                                                                                       
c  Special case:   if |indvel|=3 only atoms above zminwr get written out

      if(iabs(indout) .eq. 3) then
           nstart=1
           nend=NATMS(1)
           do itype=1,Natype
              ndumpatms(itype)=0
	      do natmchk=nstart,nend
                zchk=(RA(3*natmchk))*ax
	        if(zchk.gt.zminwr) ndumpatms(itype)=ndumpatms(itype)+1
	      enddo
              nstart=nend+1
              nend=nend+NATMS(itype+1)
	   enddo

          WRITE(lun,120) (ndumpatms(I),I=1,NATYPE)

        else
c   The usual case:
           WRITE(lun,120) (NATMS(I),I=1,NATYPE)
120        FORMAT(10i5,45x,'     # OF ATOMS OF EACH TYPE')
        endif

        WRITE(lun,130) (AMASS(I),I=1,NATYPE)
130     FORMAT(10f8.3,' MASS(i)')
c                                                                                         
        isumat=0

        DO 20 I=1,NATYPE
           WRITE(lun,156)(CPHEAD(I,K),K=1,78)
156        FORMAT(1X,80A1)
           if(iabs(indout) .ne. 2) then
              WRITE(lun,157) I
157           FORMAT('   COORDINATES OF COMPONENT #',I2)
           else
              WRITE(lun,159) I
159           FORMAT('  COORDINATES OF COMPONENT #',I2,' Without PBC!')
           endif
c                                                                                         
           do 22 k=isumat+1,isumat+natms(i)
c             Special atoms are distinguished by a tag that is 1,2,3 etc. after           
c                   the coordinates.                                                      
              if(iabs(indout) .eq. 2) then
c                      write coordinates that have not been subject to PBC:               
                 x=raNoPBC(3*K-2)*cooscale
                 y=raNoPBC(3*K-1)*cooscale
                 z=raNoPBC(3*K)*cooscale
              else
c                    write coordinates of the image inside the simulation box:            
                 x=RA(3*K-2)*cooscale
                 y=RA(3*K-1)*cooscale
                 z=RA(3*K)*cooscale
              endif
              if(indout .eq. -3 .and. z .lt. zminwr) go to 22
              WRITE(lun,140) x,y,z,itag(K),K
22            continue
140        FORMAT(3g23.15,i5,i6)
           isumat=isumat+natms(i)
20         CONTINUE
c                                                                                         
c   ------------------------------------------------------------------------              
c    Write info on atoms subject to harmonic potential:                                   
        if(nharm .gt. 0) then
           write(lun,155)
155        format('   HARMONIC CONSTRAINTS:  Center and spring ',
     +            'constant for each atom with tag 2')
           do ih=1,nharm
              x=harmcntr(3*(ih-1)+1)*cooscale
              y=harmcntr(3*(ih-1)+2)*cooscale
              z=harmcntr(3*(ih-1)+3)*cooscale
              write(lun,142) x,y,z,harmspr(ih)
              enddo
142           FORMAT(3g20.10,1x,g12.4)
        endif
c   -------------------------------------------------------------------------             
c                                                                                         
c    Write the velocities:                                                                
        if(indout .ge. 0) then
        isumat=0
        DO 21 I=1,NATYPE
           WRITE(lun,158) I
158        FORMAT('   VELOCITIES OF COMPONENT #',I2)
c                                                                                         
           do 23 k=isumat+1,isumat+natms(i)
              x=VA(3*K-2)*cooscale
              y=VA(3*K-1)*cooscale
              z=VA(3*K)*cooscale
              if(iabs(indout) .eq. 2) then
                 WRITE(lun,143) x,y,z,Potpat(k),K
              else
                 WRITE(lun,141) x,y,z,K
              endif
23            continue
           isumat=isumat+natms(i)
21         CONTINUE
        endif

141     FORMAT(3g24.15,i6)
143     FORMAT(3g24.15,2x,g24.15,i6)
c                                                                                         
c  ---------------------------------------------------------------------------            
c                                                                                         
c  Error checks:                                                                          
c   Check to see if an atom is outside the simulation cell:                               
      small=ax*0.0001
      axhalf=ax*0.5
      ayhalf=ay*0.5
      azhalf=az*0.5
      xmin=-axhalf-small
      xmax=axhalf+small
      ymin=-ayhalf-small
      ymax=ayhalf+small
      zmin=-azhalf-small
      zmax=azhalf+small
      do 150 i=1,natoms
          ind=3*(i-1)
          if(ra(ind+1) .gt. xmax .or. ra(ind+1) .lt. xmin) then
                write(lunout,9999) 1,i,ra(ind+1),xmin,xmax
9999            format(/' ERROR 999 in CONOUT: '/
     +           'coo ',i2,' of atom # ',i5,
     +          ' is out of range, Ra=',g15.5/
     +          ' while min,max =',2g15.5)
                stop
          endif
          if(ra(ind+2) .gt. ymax .or. ra(ind+2) .lt. ymin) then
                write(lunout,9999) 2,i,ra(ind+2),ymin,ymax
                stop
          endif
          if(ra(ind+3) .gt. zmax .or. ra(ind+3) .lt. zmin) then
                write(lunout,9999) 3,i,ra(ind+3),zmin,zmax
cc                stop
          endif
150       continue
c                                                                                         
c   Check to see if rigid or constrained atoms have non-zero velocity:                    
c  Permafrost atoms:                                                                      
      do 151 ip=1,naperm
          ind=3*(iperm(ip)-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 998
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 998
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 998
151       continue
c  Constrined atoms:                                                                      
c   Check to see if constrained atoms have non-zero velocity in the constrained           
c   coordinates:                                                                          
c    Constrained atoms:                                                                   
      do 153 ip=1,naconstr
          ind=3*(iconstr(ip)-1)
          if(itag(iconstr(ip)) .eq. 4) then
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 993
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 993
          endif
          if(itag(iconstr(ip)) .eq. 5) then
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 993
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 993
          endif
          if(itag(iconstr(ip)) .eq. 6) then
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 993
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 993
          endif
153       continue
cFPI:                                                                                     
c    FPI chains:                                                                          
      do 152 ip=1,nFPI

         if(itagFPI(ip) .ge. 100 .and. itagFPI(ip) .lt. 200) then
c         here endatoms should be fixed.                                                  
           ind=3*(iFPI(ip)-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 996
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 996
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 996
           ind=3*(iFPI(ip)+nimFPI(ip)-1-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 995
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 995
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 995
         endif
         if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then
c         here the centroid is constrained, so the total velocity of the                  
c         ring polymer should be zero:                                                    
           vxtot=0.0
         vytot=0.0
         vztot=0.0
         do iim=1,nimFPI(ip)
           ind=3*(iFPI(ip)+(iim-1)-1)
           vxtot=vxtot+va(ind+1)
           vytot=vytot+va(ind+2)
           vztot=vztot+va(ind+3)
           sumv=sqrt(vxtot**2+vytot**2+vztot**2)
           enddo
         if(sumv .gt. 0.1e-10) 
     +       write(lunout,9992) ip,itagFPI(ip),sumv,vxtot,vytot,vztot
9992         format(/'  WARNING 992 in CONOUT:  FPI chain ',i5,
     +        '  with tag = ',i5/
     +        '  has non-zero centroid velocity,  sumv = ',g15.5/
     +        '         vx,vy,vz = ',3g15.5)

         endif
152      continue
c                                                                                         
c     -----------------------------------------------------------------------             
c                                                                                         
      RETURN
998   continue
      write(lunout,9998) ip,iperm(ip),va(ind+1),va(ind+2),va(ind+3)
9998  format(/'  ERROR 998 in CONOUT: permafrost atom # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      return
996   continue
      write(lunout,9996) ip,iFPI(ip),va(ind+1),va(ind+2),va(ind+3)
9996  format(/'  ERROR 996 in CONOUT: first atom in FPI # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      return
995   continue
      write(lunout,9995) ip,iFPI(ip),va(ind+1),va(ind+2),va(ind+3)
9995  format(/'  ERROR 995 in CONOUT: last atom in FPI # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      return
993   continue
      write(lunout,9993) iconstr(ip),va(ind+1),va(ind+2),va(ind+3)
9993  format(/'  ERROR 993 in CONIN:  constrained atom ',i5/
     +        '         has illegal velocity,  vx,vy,vz = ',3g25.5)
      return

      END
