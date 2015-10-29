c  CONIN:    
c  version f:   add the option of having all atoms in the config in chains
c               without having chain tags on them all (thereby lifting the 
c               limitation on the number of chains)
c  version f:   add constrained FPIs (tags 500-599) for NEB calcul. on 
c               clusters.
c  EAMd92 version    (last change June 26. 1993)
c   This subroutine reads in a 3dmd configuration file.
c   This version of CONIN reads a tag for each atom after the coordinates.
c   The tag indicates what kind of dynamics apply to the atom.                
c      tag #:     dynamics:
c        0        plain Newton with no constraints except PBC
c        1        permafrost atom, force and velocity are zeroed
c        2        atom subject to external harmonic potential (not fully 
c                 implemented)
c        3            
c        4        atom constrained in y and z
c        5        atom constrained in x and z    
c        6        atom constrained in x and y    
c        7        atom constrained in x
c        8        atom constrained in y
c        9        atom constrained in z
c       10        atom subject to random force (Langevin dynamics)
c       11        H2O molecule with center of mass fixed.
c       12        H2O molecule with center of mass fixed in x direction.
c       13        H2O molecule with center of mass fixed in y direction.
c       14        H2O molecule with center of mass fixed in z direction.
c       15        H2O molecule with CM fixed on the vertical plane at 30degree
c       16        H2O molecule with CM fixed on the vertical plane at 150degree
c       17        H2O molecule with CM fixed in y and z (just moves along x)
c       18        H2O molecule with CM fixed in x and z (just moves along y)
c       19        H2O molecule with CM fixed in x and y (just moves along z)
c     100-199     open FPI chain with fixed endpoints (useful for barrier 
c                 calculations)
c     200-299     circular FPI (the `ordinary' diagonal density matrix 
c                 element FPI)
c     300-399     circular FPI with fixed centroid
c     500-599     FPI with constraints (second digit analogous to above)
c         500-539       all images fixed (permafrost chain)
c         540-549       images constrained in y and z
c         550-559       images constrained in x and z    
c         560-569       images constrained in x and y    
c         570-579       images constrained in x
c         580-589       images constrained in y
c         590-599       images constrained in z
cAllinCh:
c     100x        If the first N atoms in the configuration file have
c                 a tag of 1000, then all atoms in the configuration file
c                 will be interpreted as images in chains of length N.
c                 If the first tag is 1001 then the first chain is all 
c                 permafrost.
c                 If the first tag is 1004 then the images of the first chain
c                       are constrained in y and z,   etc.
c                 Subsequent chains can have tags of 0, 1, 2, 3, etc.
c                 but all images in a given chain must have the same tag.
c                 Internally in the program, the tag on atoms in the first
c                 chain is set to x, and then changed to 100x before writing 
c                 co.con
c     110x        Same as 100x, except that that the chains will have fixed
c                 endatoms, the tags of all endatoms are set to 1.
c-----------------------------------------------------------------------------


      SUBROUTINE CONIN
c 
      implicit real*8 (a-h,o-z)

c
c             ax, ay and az are the lengths of the simulation cell sides.
c             alpha is the angle between the ax and ay vectors, beta the
c             angle between the ay and az and gamma betwee ax and az.

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comluns.cmn'

      COMMON /ISEEDS/ ISEED1

      data rdprdg / 0.1745329252d-01 /


c-----------------------------------------------------------------------

      indallinchains=0
c          by default, atoms are not all images in chains.
      indallfixend=0

      READ(lunci,*,end=91,err=91) ISEED1
c                  Seed for the random number generator. 
      if(iseed1 .eq. 0) then
            write(6,*) '    WARNING from CONIN:   iseed = 0 '
            write(lunout,*) '    WARNING from CONIN:   iseed = 0 '
      endif
      go to 101
91    continue
      write(lunout,*)' CONIN: end or err in reading ISEED (empty file?)'
      write(6,*) '   CONIN:  end or err in reading ISEED, (empty file?)'
      return

101   continue

100   FORMAT(I20)

      read(lunci,*) time0
c               this is the simulated time initially (simply for reference).
              
      read(lunci,*) ax,ay,az

        ndimens=3
        if(az .lt. 0.1e-07) then
c         since the az dimension of the box is nearly zero it is assumed the
              
c         simulation is two dimensional.  Volume then actually means area:
                
          ndimens=2
c         its safer to make az finite but small since
                                     
c                 PBC will still be applied in z:
                                         
          az=1.0e-10
        endif
        volm0=ax*ay*az
        volm=volm0
        ax0=ax
        ay0=ay
        az0=az
c                store the initial values of the box sidelengths and volume.
              

      read(lunci,*) alpha,beta,gamma
c               ax, ay and az are the sides of the simulation cell,
                       
c               alpha the angle between ax and ay, beta between ay and az, etc
c DEVELOPMNT:
                                                                             
c     At this stage only 90 deg angles are allowed.
                                       
        diff=abs(alpha-90)+abs(beta-90)+abs(gamma-90)
        if(diff .gt. 0.00001) then
           write(lunout,*) ' WARNING: all box angles must be 90 for now'
           alpha=90.0d+00
           beta=90.0d+00
           gamma=90.0d+00
        endif

      read(lunci,*) pmass,volv
c         pmass   is the piston mass in isobaric simulations.
                             
c         volv    is the time derivative (generalized velocity) of the volume.
            

      read(lunci,*) naperm,nharm
c          naperm is the number of rigid (permafrost) atoms.  This value is
               
c                   compared with the number of atoms tagged with 1 and the
               
c                   program stops if there is an error.
                                   
c          nharm  is the number of atoms that are subject to external harmonic
            
c                   potential. This value is
                                              
c                   compared with the number of atoms tagged with 2 and the
               
c                   program stops if there is an error. 
                                  

      sinalp=sin(rdprdg*alpha)
      cosalp=sqrt(1.-sinalp**2)
      sinbet=sin(rdprdg*beta)
      cosbet=sqrt(1.-sinbet**2)
      sinbet=sin(rdprdg*beta)
      cosbet=sqrt(1.-sinbet**2)
c DEVELOPMNT:
                                                                             
c   The following three lines are for 90deg angles only
                                   
c   and need to be generalized:
                                                           
      axhalf=ax*0.5
      ayhalf=ay*0.5
      azhalf=az*0.5

c  READ number of atoms, mass of atoms, coordinates and velocities
                        
c  for each type of atoms:

      READ(lunci,*) NATYPE
c                 number of types of atoms (components).

      IF(NATYPE.GT.0) THEN
        READ(lunci,*) (NATMS(I),I=1,NATYPE)
        READ(lunci,*) (AMASS(I),I=1,NATYPE)
        do ic=1,natype
         if(amass(ic) .le. 0.1e-30) then
           write(lunout,*) '  ERROR from CONIN:  mass is zero, ic = ',ic
         endif
      enddo
c 

        NATOMS=0
        DO 10 I=1,NATYPE
10         NATOMS=NATOMS+NATMS(I)

        indst=0
cFPI: 
        nsumFPI=0
        itagFPIcur=-999
        do i=100,MAXFPI
          numFPI(i)=0
          nimFPI(i)=0
        enddo
        do i=1,10
          nFPIcom(i)=0
        enddo
 
        DO 20 I=1,NATYPE
            Read(lunci,156) (CPHEAD(I,K),K=1,80)
156         FORMAT(80A1)
            Read(lunci,100)
c                  skip line to allow for header saying COORDINATES ... 
                  
            do 21 ka=1,natms(i)
cw            write(6,*) '  conin:  ka=',ka 
                                                                  
                ind=3*(ka-1+indst)
                read(lunci,*) x,y,z,itag(ka+indst)
cw                write(6,*) '  CONIN: ka=',ka,' x,y,z = ',x,y,z,
                         
cw     +                 '  tag = ',itag(ka+indst)
                                        
c             Check tags:                                      
cAllinCh:  ---------
              if(I .eq. 1 .and. ka .eq. 1) then
cw                 write(6,*) ' check out the tag on first atom'
                 ifirsttag=itag(1)
cw                 write(6,*)'  ifirsttag,itag(1)=',ifirsttag,itag(1)
                 if(itag(1).ge.1000.and.itag(1).le.1200) then
c               all atoms in this configuration will be assumed images 
c               in chains
                   indallinchains=1
                   ifindlength=1
                   iFPI(1)=1
                   if(itag(1).ge.1100) indallfixend=1
                 endif
              endif
c                      close `if(I .eq. 1 .and. ka .eq. 1)'

              if(ka .gt. 1 .and. ifindlength .eq. 1) then
c           find how long the first chain is (and this becomes the 
c           length of all ch)
                 if(itag(ka+indst) .ne. ifirsttag) then
cw                    write(6,*)' Find beginning of second chain, ka=',ka
cw                    write(6,*)' itag,ifirsttag=',itag(ka+indst),ifirsttag
                    ifindlength=0
c                      mark that an atom with tag different from the tag of
c                      the first atom has been found.
                    nimFPI(1)=ka-1
                    if(ka .le. 1) then
                      write(6,*) '  ERROR in conin:  nimFPI(1)=0, ka=',ka
                      stop
                    endif
                    nsumFPI=natoms/nimFPI(1)
                    nFPI=nsumFPI
                    if(indallfixend .eq. 1) then
                         itag(nimFPI(1))=1
                         itag(1)=1
                    endif
                    do icomp=1,10
                       nFPIcom(icomp)=Natms(icomp)/nimFPI(1)
                    enddo
                    do ifa=1,nFPI
                      iFPI(ifa)=1+nimFPI(1)*(ifa-1)
                      nimFPI(ifa)=nimFPI(1)
                    enddo
                    if(natoms .ne. nFPI*nimFPI(1)) then
                      write(6,*) '  ERROR in conin: allinchains, but'
                      write(6,*) '  nFPI,nimFPI(1)=',nFPI,nimFPI(1)
                      stop
                    endif
                 endif
              endif

              if(indallinchains .eq. 1 .and. ifindlength .eq. 0) then
                n=(ka+indst-1)/nimFPI(1)
                idiff=ka+indst-nimFPI(1)*n
                if(idiff .eq. 2) itagFPI(n+1)=itag(ka+indst)
                if(indallfixend .eq. 1) then
c              set the tag of endpoint atoms to 1:
                 if(idiff .eq. 1 .or. idiff .eq. nimFPI(1)) then
                    itag(ka+indst)=1
                 end if
                endif
              endif
              if(indallfixend .eq. 1 .and. itag(ka+indst) .ge. 1100)
     +              itag(ka+indst)=itag(ka+indst)-1100
cERB
              if(indallfixend .eq. 1 .and. ka+indst .eq. 1)
     +              itag(ka+indst)=1
              if(indallinchains .eq. 1 .and. itag(ka+indst) .ge. 1000)
     +              itag(ka+indst)=itag(ka+indst)-1000

cw              write(6,*)' endchecks on AllinCh,itag=',itag(ka+indst)


c          Stick the coordinates of all the atoms included into array RA:

                if(itag(ka+indst) .ge. 0) then
                     ra(ind+1)=x
                     ra(ind+2)=y
                     ra(ind+3)=z
                        aveRA(ind+1)=ra(ind+1)
                        aveRA(ind+2)=ra(ind+2)
                        aveRA(ind+3)=ra(ind+3)
c             the array aveRA is used to calculate the average coordinates.

                endif
                if(itag(ka+indst) .eq. -1) then
c                     atoms tagged with -1 are not included in the calculation            
c                     but the coordinates are stored in the array 'rinfraperm'.           
                     if((ninfraperm+1) .gt. MAXCOO) then
                      write(lunout,*) 'ERROR CONIN: ninfraperm>MAXCOO'
                        stop
                     endif
                     rinfraperm(3*ninfraperm+1)=x
                     rinfraperm(3*ninfraperm+2)=y
                     rinfraperm(3*ninfraperm+3)=z
                     ninfraperm=ninfraperm+1
                     if(I .eq. 1) nAinfraperm=nAinfraperm+1
                endif


21              continue
            indst=natms(i)+indst
20          continue

cendAllinCh.  -------------------------

            indst = 0
            nsumperm=0
            nsumharm=0
            ninfraperm=0
            nAinfraperm=0
            nsumconstr=0
            nlangevin=0

        DO 40 I=1,NATYPE
            do 41 ka=1,natms(i)
                ind=3*(ka-1+indst)


c          Check if permafrost:
                if(itag(ka+indst).eq.1 .or.
     +            (itag(ka+indst).ge.500.and.itag(ka+indst).le.539))then
c                   this is a permafrost atom or image:
                                            
                    nsumperm=nsumperm+1
                    iperm(nsumperm)=ka+indst
                endif
                if(itag(ka+indst) .eq. 2) then
c                   this is an atom subject to external harmonic potential:
               
                    nsumharm=nsumharm+1
                    iharm(nsumharm)=ka+indst
                endif
                if((itag(ka+indst) .eq. 4) .or.
     +          (itag(ka+indst).ge.540.and.itag(ka+indst).le.549))then
c                   this is an atom constrained in y and z:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if((itag(ka+indst) .eq. 5) .or.
     +          (itag(ka+indst).ge.550.and.itag(ka+indst).le.559))then
c                   this is an atom constrained in x and z:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if((itag(ka+indst) .eq. 6)  .or.
     +          (itag(ka+indst).ge.560.and.itag(ka+indst).le.569))then
c                   this is an atom constrained in x and y:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if((itag(ka+indst) .eq. 7) .or.
     +          (itag(ka+indst).ge.570.and.itag(ka+indst).le.579))then
c                   this is an atom constrained in x:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if((itag(ka+indst) .eq. 8)  .or.
     +          (itag(ka+indst).ge.580.and.itag(ka+indst).le.589))then
c                   this is an atom constrained iny:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if((itag(ka+indst) .eq. 9)  .or.
     +          (itag(ka+indst).ge.590.and.itag(ka+indst).le.599))then
c                   this is an atom constrained in z: 
                              
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if(itag(ka+indst) .eq. 10) then
c                   this is an atom subject to random force (Langevin 
c                   dynamics):                               
                    nlangevin=nlangevin+1
                    ilangevin(nlangevin)=ka+indst
                endif

cERB Constrains for the CM of an H2O molecule added by ERB
                if((itag(ka+indst) .ge. 11) .or. 
     $               (itag(ka+indst) .le. 19))  then
c                   this is oxygen indicates a molecule with the center
c                   of mass fixed  
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
cFPI:   
              if(itag(ka+indst).ge.100 .and. itag(ka+indst).lt.600)then
c                 this is an atom in a FPI chain:
                                         
c                 NB! it is assumed here that all images corresponding to 
                
c                 a given FPI chain are consecutive in ci.con.
                            
                if(itag(ka+indst).ge.500.and.itag(ka+indst).lt.600)then 
c                   this is an image constrained in some way:
                               
                    nsumconstr=nsumconstr+1
                    iconstr(nsumconstr)=ka+indst
                endif
                if(itag(ka+indst) .ne. itagFPIcur) then
c                     start a new FPI chain: 
                                             
c                       first, check whether this tag has been seen before:
               
                   do ich=1,nsumFPI
                    if(itagFPI(ich) .eq. itag(ka+indst)) go to 994
                   enddo
                      itagFPIcur=itag(ka+indst)
                      nsumFPI=nsumFPI+1
                   numFPI(itag(ka+indst))=nsumFPI
                      itagFPI(nsumFPI)=itag(ka+indst)
                      iFPI(nsumFPI)=ka+indst
                      nFPIcom(i)=nFPIcom(i)+1
                   nimFPI(nsumFPI)=1
                else
c                       this is a new image in the current FPI chain:
                     
                  nimFPI(nsumFPI)=nimFPI(nsumFPI)+1
                endif
              endif


              numIMAGE=nimFPI(1)

c-----------------------------------------------------------------------
41              continue
            indst=natms(i)+indst
40          continue

c
         if(naperm .ne. nsumperm) then
            write(lunout,8996) naperm,nsumperm
8996     format(/'# WARNING from CONIN: naperm .ne. nsumperm,  naperm ='
     +                 ,i5,' nsumperm=',i5)
            naperm=nsumperm
         endif
         if(nharm .ne. nsumharm) then
            write(lunout,8995) nharm,nsumharm
8995        format(/'  WARNING from CONIN:   nharm .ne. nsumharm, ',
     +                 'nharm =',i5,' nsumharm=',i5)
            nharm=nsumharm
         endif

         naconstr=nsumconstr

cFPI:
         nFPI=nsumFPI
         if(indallfixend .eq. 1) then
           do ifa=1,nFPI
            itag(1+(ifa-1)*nimFPI(1))=1
            itag(ifa*nimFPI(1))=1
           enddo
         endif

c-----------------------------------------------------------------------
c   Read info on atoms subject to external harmonic potential:
                            
      if(nharm .gt. 0) then
       read(lunci,*)
c               skip header 
                                                              
       do ih=1,nharm
          ind=3*(ih-1)
          read(lunci,*,end=98,err=98)
     +      harmcntr(ind+1),harmcntr(ind+2),harmcntr(ind+3),harmspr(ih)
       enddo
       go to 96
98     continue
       write(lunout,*)  'ERROR in CONIN:  end or err in reading harmonic
     + info'
       stop
      endif
c             close 'if(nharm .gt. 0)'
                                                    

c-----------------------------------------------------------------------

96      continue

c   Read the velocities:
                                                                  
         indst=0
         do 23 i=1,natype
c                         NB! CPHEAD is not read for velocities. 
                         
           Read(lunci,100,end=99,err=99)
c                 skip line to allow for header saying VELOCITIES ... 
                    
c       First read one line to test the waters. Check whether the 
c       potential per atom is given:
                                                                    
cc  ***  This doesnt seem to work any more ***
cc           indpot=1
cc           vax=0.
cc           vay=0.
cc           vaz=0.
cc           read(lunci,*,end=97,err=97) vax,vay,vaz,poten
cc             va(3*indst+1)=vax
cc             va(3*indst+2)=vay
cc             va(3*indst+3)=vaz
c            check whether the value of poten is simply the atom number:                  
cc             if(abs(poten-indst-1.) .lt. 0.1e-06) then
cc               poten=0.0
cc               indpot=0
cc             endif
cc             Potpat(indst+1)=poten

cc             write(6,*) '   poten(',indst+1,') = ',poten                                

cc             go to 102
cc97         continue
c         if END or ERR comes up in the first read statement,                             
c         then assume the fourth column (Potpat) is missing:                              

cc             write(6,*)'   CONIN: potential per atom not found, va =',
cc     +                  vax,vay,vaz,'  indst = ',indst,'  itype= ',i

cc             indpot=0
cc             va(3*indst+1)=vax
cc             va(3*indst+2)=vay
cc             va(3*indst+3)=vaz
cc102        continue

cc           write(6,*) 'indpot = ',indpot,'     for component ',i                        

cc           do 22 ka=2,natms(i)

           indpot=0
           do 22 ka=1,natms(i)
             ind=3*(ka-1+indst)
             indat=ka+indst
             if(indpot .eq. 1) then
                        read(lunci,*,end=95,err=95)
     +                   va(ind+1),va(ind+2),va(ind+3),Potpat(indat)
             else
                               read(lunci,*,end=95,err=95)
     +                              va(ind+1),va(ind+2),va(ind+3)
             endif
22           continue
           go to 197
c                                                                                         
99         continue
c    if END or ERR in reading velocities, then assume velocities are all zero:            
           write(6,220) i
220        format('    CONIN: velocity header not found, set va=0. for',
     +               ' all atoms of component',i5)
            go to 104
95          continue
            write(6,*)'    CONIN: velocities not found, set',
     +                ' velocity = 0. for all atoms of component',i
104         continue
            do 198 ka=1,natms(i)
               if(indpot .eq. 0) Potpat(ka+indst)=0.0
               ind=3*(ka-1+indst)
               va(ind+1)=0.0
               va(ind+2)=0.0
               va(ind+3)=0.0
198            continue
c 

 197           continue
            indst=natms(i)+indst
 23      CONTINUE
      ENDIF
c 
c-----------------------------------------------------------------------
c
      RETURN

994   continue
      write(lunout,9994) nsumFPI,itagFPI(ich),itag(ka+indst)
9994  format(/'  ERROR 994 in CONIN:  FPI images not adjacent ',
     +        '  nsumFPI = ',i4,' itagFPI(ich) = ',i4,
     +        '  itag(ka+indst) = ',i4)
      stop

      END















