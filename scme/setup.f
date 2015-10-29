c     THIS SUBROUTINE PREPARES THE INPUT CONFIGURATION FOR DYNAMICS                       
c          - Characterizes and checks the interaction potential
c          - Does basic checks on the input configuration
c          - Scales the coordinates and velocities.                                       
c          - CALLS RANVEL - WHICH ASSIGNS RANDOM VELOCITIES                               
c          - CALLS FORCE - WHICH INITIALIZES F ARRAY FOR DYNAMICS                         
c          - COMPUTES INITIAL (REFERENCE) VALUES                                          
c          - CHARACTERIZES THE STARTING CONFIGURATION                                     
c                                                                                         
      SUBROUTINE SETUP(eQM,dEQM,FAOUT,EPOTOUT,ETOUT,natm,QPOLEOUT)
c                                                                                         
      implicit real*8 (a-h,o-z)
      integer natm
      real*8 eQM(3,natm)
      intent(in) eQM
      real*8 dEQM(3,3,natm)
      intent(in) dEQM
      dimension FAOUT(*), ETOUT(*), QPOLEOUT(*)
      intent(out) FAOUT, ETOUT, QPOLEOUT
      real*8 EPOTOUT
      intent(out) EPOTOUT 

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comtgr.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comoutcntrl.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/compotent.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comdeposit.cmn'
      include '../commonblks/comintlis.cmn'
      include '../commonblks/constraints.cmn'

      COMMON /ISEEDS/ ISEED1
      COMMON /REFENC/ REFPRE,REFVOL,REFPOT,REFKIN,
     X                REFENT,REFENG,temp0,refpis

c      The following commonblock is here just for testing the potential:
cc      common / neighborlist /
cc     +               r2pr(150000),delpr(500000),r2st(2002),
cc     +               rhoij(150000),rhofij(500000),rhofji(500000),
cc     +               rhoji(150000),potpr(250000),
cc     +               drhoijr(150000),drhojir(150000),
cc    +               indpra(2,150000),iptpr1(2002,3),iptcom(4),
cc     +               nintst(3),natm1(3),natm2(3),nintp(3)

      dimension ra1pl(3),ra2pl(3),fa1pl(3),
     +          fa2pl(3),virpl(6),eldensat(2),natmsst(MAXCOMP),
     +          rpott(MAXPOTCH),pottes(MAXPOTCH),FORTES(MAXPOTCH),
     +          rhotes(MAXPOTCH),fofrhotes(MAXPOTCH),enetes(MAXPOTCH)
      dimension embforFint(MAXATOMS,3),natm1forFint(3),
     +          rho1forFint(MAXATOMS,3)

      logical Iskippottbl
      character*80 dummyCH
c-----------------------------------------------------------------------
cERB Added for general constraint
      if (1 .lt. 0) then
      
      if (genConstraints) then
         open(99,file = 'constraints.dat')
c$$$         read(99,'(A)') dummyCH
c$$$         read(99,'(A)') dummyCH
         read(99, *) bConst
         A2 = 0.d0
         do i = 1, natoms
            read (99,*) (constraint(3*(i-1)+j), j=1,3), nn         
            A2 = A2 + constraint(3*(i-1)+1) * constraint(3*(i-1)+1)
            A2 = A2 + constraint(3*(i-1)+2) * constraint(3*(i-1)+2)
            A2 = A2 + constraint(3*(i-1)+3) * constraint(3*(i-1)+3)
         end do
         if (A2 .ge. 1.d-15) then
            A2 = sqrt(A2)
            do i = 1, natoms
               do j = 1, 3
                  constraint(3*(i-1)+j) = constraint(3*(i-1)+j) / A2
               end do
            end do
         end if
         close(99)
      end if

c      end if

      if (1 .lt. 0) then

      if (genConstraints) then
         open(99,file = 'constraints.dat')
         read(99, *) ibConst
         A2 = 0.d0
         do i = 1, nAtms(2)
            read (99,*) (constraint(3*(i-1)+j), j=1,3), nn         
            A2 = A2 + constraint(3*(i-1)+1) * constraint(3*(i-1)+1)
            A2 = A2 + constraint(3*(i-1)+2) * constraint(3*(i-1)+2)
            A2 = A2 + constraint(3*(i-1)+3) * constraint(3*(i-1)+3)
         end do
         if (A2 .ge. 1.d-15) then
            A2 = sqrt(A2)
            do i = 1, natoms
               do j = 1, 3
                  constraint(3*(i-1)+j) = constraint(3*(i-1)+j) / A2
               end do
            end do
         end if
         close(99)
      end if
      
      end if
      end if
cERB end of constraints
c-----------------------------------------------------------------------


      Iskippottbl=.true.

      maxnpr=0
      maxneb=0

                                                                                         
c   Find the 'skin depth', i.e. distance between potential cutoff and the                 
c   cutoff for the list of interacting pairs.                                             
        skndpt(1)=rskin(1)-rcut(1)
        skndpt(2)=rskin(2)-rcut(2)
        skndpt(3)=rskin(3)-rcut(3)

        rskinmax=dmax1(rskin(1),rskin(2),rskin(3))

c      -------------------------------------------------------------------                

c  Interaction Potential:                                                                 

      npotint=natype*(natype+1)/2
      if(natype .eq. 1) write(lunout,369)
369   format('#'/
     +     '# --------------------------------------------------------'/
     +     '# INTERACTION POTENTIAL AND FORCE:'/
     +     '#',14x,'                          A-A:  ')
      if(natype .eq. 2) write(lunout,366)
366   format('# ------------------------------------------------------'/
     +       '# INTERACTION POTENTIAL AND FORCE:'/'#',
     +     14x,'                          A-A:        B-B:       A-B:')
      write(lunout,362) (rcut(k),k=1,npotint)
362   format('#          cutoff distance, rcut =',3f10.4)
      write(lunout,368) (skndpt(k),k=1,npotint)
368   format('#          neighbor list, skndpt =',3f10.4)
      write(lunout,413)
413   format('#                       Parameters:')

      do ip=1,npotpar
         write(lunout,365)
cc     +                      (potheader(jk,ip),jk=1,30),                                 
     +                      (potpar(ip,npk),npk=1,npotint)
      enddo
365   format('#              ',
cc     +        30a1,' = ',                                                               
     +                        19x,3f12.4)
      write(lunout,412)
412   format('#')

      mxipot=npotint                                                                    

c ------------------------------------------------------------------------                

c  CHECK THE POTENTIAL AT A FEW POINTS and PRINT OUT A TABLE OF VALUES:                   
      if(Iskippottbl) go to 911
      indupd=1
c   atom 2 is at the center, while atom 1 is displaced along x-axis.                      
      ra1pl(2)=0.
      ra1pl(3)=0.
      ra2pl(1)=0.
      ra2pl(2)=0.
      ra2pl(3)=0.

c   store various info that gets scrambled in mimicking dimer:
      natomsst=natoms
      natoms=2
      do i=1,10
        natmsst(i)=natms(i)
      natms(i)=0
      enddo
      itagst1=itag(1)
      itagst2=itag(2)
      itag(1)=0
      itag(2)=0
      numIMAGEst=numIMAGE
      numIMAGE=1


      WRITE(lunout,361)                                             
361   format('#'/'#'/
     +       '#   Evaluate the potential for dimer at a few points:'/
     +       '#   -------------------------------------------------')

      DO 376 i=1,npotint                                                                
c          loop over potential types:                                                     
        indupd=1                                                                        

c    Find how many points, spaced by 0.5, are needed to cover the intereaction
c    range:
        npotch=rskinmax*2.0+0.51
        if(npotch .gt. maxpotch) npotch=maxpotch

        DO 382 j=1,npotch                                                                    
          ra1pl(1)=j*0.5                                                                
          rpott(j)=ra1pl(1)                                                             

cPres:      In EAM version the coordinates should NOT be scaled by ax                     
c           before calling GAGAFE.                                                        

          do if=1,6                                                                     
	    fa(if)=0.0                                                                         
	  enddo                                                                                
          rhoij(1)=0.0
          rhoji(1)=0.0

          if(i .eq. 1) then
             it1=1
             it2=1
             n1=2
             n2=2
             isame=1
          endif
          if(i .eq. 2) then
             it1=2
             it2=2
             n1=2
             n2=2
             isame=1
          endif
          if(i .eq. 3) then
             it1=1
             it2=2
             n1=1
             n2=1
             isame=0
          endif
          iargFA2=4-3*isame
cERB          CALL GAGAFE(n1,ra1pl(1),FA(1),n2,ra2pl(1),FA(iargFA2)
cERB     +                ,utot,virpl,i,isame)                                                  
          POTTES(J)=utot                                                                

cPres:    get actual forces from GAGAFE in the EAM version and                            
c         no need to rescale.                                                    

          FORTES(J)=fa(1)                                                               

	  RHOTES(j)=rhoij(1)          
          rho1forFint(1,it1)=RHOTES(j)
          natm1forFint(it1)=1
          call FINT(rho1forFint,natm1forFint,it1,embforFint)
          fofrhotes(j)=embforFint(1,it1)

	  RHOTES(j)=rhoji(1)
          rho1forFint(1,it2)=RHOTES(j)
          natm1forFint(it2)=1
          call FINT(rho1forFint,natm1forFint,it2,embforFint)
          fofrhotes(j)=fofrhotes(j)+embforFint(1,it2)

          enetes(j)=fofrhotes(j)+POTTES(j)                                                 

382       CONTINUE                                                                      

cx          CALL EMBED(FRHOTOT,embvir)                                                    
cx          fofrhotes(j)=0.5*FRHOTOT                                                      
cx	    RHOTES(j)=eldensat(1)                                                                
cx	    enetes(j)=2.0*fofrhotes(j)+POTTES(j)                                                 


c       write results for dimers:                                                           
          WRITE(lunout,363) i
363       format('#'/'# Potential ',i3,':'/
     +       '#  Distance:   Pair Pot:   El.Dens:    Embed.En:',
     +       ' Force (phi):  Tot.Energy:')                                                         

          do jj=1,npotch
             WRITE(lunout,370) Rpott(JJ),POTTES(JJ),RHOTES(JJ),
     +                        FofRHOTES(JJ),FORTES(JJ),enetes(JJ)                                    
370          format('#p ',f4.1,'  ',5g13.4)
          enddo
          write(lunout,410) 
410       format('#p')

376     CONTINUE                                                                        


c    Print cutoff values and check whether potential and electron density go
c    to zero beond the cutoff.

       WRITE(lunout,364) 
364    format('#'/'#    Cutoff values:'/
     +        '#   Potential:      Phicut:      rho1cut:      rho2cut:')
       do ip=1,npotint
          WRITE(lunout,367) ip,phicutst(ip),rhocut1st(ip),rhocut2st(ip)
367       format('#   ',i3,'     ',3g15.5)

          if(Rpott(npotch) .gt. rcut(ip)) then
             if(abs(POTTES(15)) .gt. 0.1e-12) then
               write(6,*) '  ERROR:  phi doesnt go properly to zero'
               write(lunout,*)' ERROR:  phi doesnt go properly to zero'
             endif
             if(abs(RHOTES(15)) .gt. 0.1e-12) then
               write(6,*)' ERROR:  rho doesnt go properly to zero'
               write(lunout,*)' ERROR: rho doesnt go properly to zero'
             endif
          endif
       enddo


c                         -------------

c    Restore the tags:                                                                    
      itag(1)=itagst1
      itag(2)=itagst2
      natoms=natomsst
      do i=1,10
        natms(i)=natmsst(i)
      enddo
      numIMAGE=numIMAGEst

911     continue

c  --------------------------------------------------------------------------             
c  --------------------------------------------------------------------------             
c                                                                                         
c   CHECK THE PROPERTIES OF THE STARTING CONFIGURATION and WRITE INFO:                          
c                                                                                         

c   If the tag is 0 there is no special treatment.                                        
c    -  -   -  -  1 the atom is rigid (a permafrost atom).                                
c    -  -   -  -  2 the atom subject to external harmonic potential.                      
c    -  -   -  -  3 the atom is not included in the calculation.                          
c    -  -   -  -  4 the atom is constrained, can only move in x direction                 
c    -  -   -  -  5 the atom is constrained, can only move in y direction                 
c    -  -   -  -  6 the atom is constrained, can only move in z direction                 
c    -  -   -  -  100 to 199 the atom is an image in FPI chain.                           
c                          The first and last atoms with each of these tag                
c                          values, say 100, is treated as rigid.                          
c    -  -   -  -  200 to 299 the atom is an image in FPI chain loop.                      
c                          The first and last atoms with each of these tag                
c                          values, say 10, is treated as being connected.                 
c    -  -   -  -  300 to 399 the atom is an image in FPI chain loop.                      
c                          The first and last atoms with each of these tag                
c                          values, say 10, is treated as being connected.                 
c                          The centroind is constrained to be in the initial              
c                          position by subtracting the total force from the               
c                          force on each image and subtracting the total velocity         
c                          from the velocity of each image.                               
c                                                                                         
c                                                                                         
c             The coordinates should be Cartesian and range from                          
c                x = [-0.5*ax , 0.5*ax], etc.                                             
c             i.e. the x-axis is parallel to the ax side of the cell.                     
c             The unit of length is determined by the interaction potential               
c             (for L-J it is usually SIGMA).                                              
c             The time (in calculating velocities) should be in units of                  
c             TAU (= [length]*sqrt([mass]/[energy]).                                      
c             The unit of energy is determined  by the interaction potential              
c             (for L-J it is usually SIGMA) and the unit of mass is derived               
c             from the mass given in the input configuration.                             

       axhalf=0.5*ax
       ayhalf=0.5*ay
       azhalf=0.5*az

c      print *,ax
      write(lunout,308) ax,ay,az,alpha,beta,gamma
308   format('#'/
     + '# ---------------------------------------------------------'/
     +       '# Simulation cell:               ax = ',g12.4/
     +        '#                                ay = ',g12.4/
     +        '#                                az = ',g12.4/
     +        '#                             alpha = ',g12.4/
     +        '#                             beta  = ',g12.4/
     +        '#                             gamma = ',g12.4)
c                                                                                         

      if(natoms .le. 0) then
         write(6,*) '   SETUP:    No atoms, natoms = ',natoms
         write(lunout,*) '   SETUP:    No atoms, natoms = ',natoms
         stop
      endif
      IF(NATOMS.GT.maxatoms) THEN
          WRITE(lunout,991) natoms
991       format(//'  **** ERROR * *** ',
     +          'TOO MANY ATOMS.  NATOMS =',I10)
          STOP
      ENDIF

      if(natype .gt. 2) then
       write(lunout,*)' ERROR SETUP: natype=',natype,'  max 2 allowed'
       stop
      endif
      WRITE(lunout,310) NATYPE
310   FORMAT('#'/'#  Number of types of atoms (components): ',i3/
     +        '#      Type:   Number of atoms:  Mass of each atom:')
      WRITE(lunout,311) (I,NATMS(I),AMASS(I),I=1,NATYPE)
311   format('#     ',i5,'            ',i7,10x,g12.4)

      if(IPRES) WRITE(lunout,312) pmass
312   format('#'/'#  Piston Mass = ',g15.4)
      write(lunout,412) 

cFPI:                                                                                     
      if(nFPI .gt. 0) then
        WRITE(lunout,314) nFPI,nFPIcom(1),nFPIcom(2)
314     format('#'/'#  Total of ',i4,' Feynman Path Integrals,',i5,
     +        ' of component 1,',i3,' of component 2'/'#',
     +        ' # of FPI:     tag:     type:    first im: ',
     +        ' #of images:  spring cons:',
     +        '  numFPI:')
       nstart=0
       nfinish=0
       do itype=1,natype
        nstart=nfinish+1
        nfinish=nstart+nFPIcom(itype)-1
        do if=nstart,nfinish
cAllinCh:
          if(indallinchains .eq. 0) then
           write(lunout,315) if,itagFPI(if),itype,iFPI(if),nimFPI(if),
     +                       sprconFPI(itype),numFPI(itagFPI(if))
          else
           write(lunout,315) if,itagFPI(if),itype,iFPI(if),nimFPI(if),
     +                       sprconFPI(itype)
          endif
315        format('# ',i8,i9,i6,i13,i13,'   ',g12.4,i9)
        enddo
       enddo
       write(lunout,412)
      endif


c  Do various error checks:   (Move this from conin, july '92)


c      print *, itag

c   Check to see if an atom is outside the simulation cell:
      small=ax*0.0001
      xmin=-axhalf-small
      xmax=axhalf+small
      ymin=-ayhalf-small
      ymax=ayhalf+small
      zmin=-azhalf-small
      zmax=azhalf+small
      do 150 i=1,natoms
          ind=3*(i-1)
          if(itag(i) .lt. 0 .or. itag(i) .gt. 600) go to 990
          if(ra(ind+1) .gt. xmax .or. ra(ind+1) .lt. xmin) then
                write(lunout,9999) 1,i,ra(ind+1),xmin,xmax
9999            format(/' ERROR 999 in SETUP: coo ',i2,' of atom # ',i5,
     +          ' is out of range, Ra=',g12.4,' while min,max =',2g12.4)
                stop
          endif
          if(ra(ind+2) .gt. ymax .or. ra(ind+2) .lt. ymin) then
                write(lunout,9999) 2,i,ra(ind+2),ymin,ymax
                stop
          endif
          if(ra(ind+3) .gt. zmax .or. ra(ind+3) .lt. zmin) then
                write(lunout,9999) 3,i,ra(ind+3),zmin,zmax
                stop
          endif
150       continue

c   Check to see if rigid atoms have non-zero velocity:
c    Permafrost atoms:                                                                    
      do 151 ip=1,naperm
          ind=3*(iperm(ip)-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 996
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 996
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 996
151       continue
c                                                                                         
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
      nFPIstr=0
      nFPIcir=0
      nFPIcnt=0
c         count the various types of FPI chains.

      do 152 ip=1,nFPI

         if((itagFPI(ip) .ge. 100 .and. itagFPI(ip) .lt. 200) .or.
     +      (itagFPI(ip) .ge. 500 .and. itagFPI(ip) .lt. 600)) then
c         here endatoms should be fixed.                                                  
           nFPIstr=nFPIstr+1
           ind=3*(iFPI(ip)-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 989
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 989
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 989
           ind=3*(iFPI(ip)+nimFPI(ip)-1-1)
            if(abs(va(ind+1)) .gt. 0.1e-20)  go to 995
            if(abs(va(ind+2)) .gt. 0.1e-20)  go to 995
            if(abs(va(ind+3)) .gt. 0.1e-20)  go to 995
         endif
         if(itagFPI(ip) .ge. 200 .and. itagFPI(ip) .lt. 300) then
c          here the FPI is a circular one that is free to move.
           nFPIcir=nFPIcir+1
         endif
         if(itagFPI(ip) .ge. 300 .and. itagFPI(ip) .lt. 400) then
c         here the centroid is constrained, so the total velocity of the                  
c         ring polymer should be zero:                                                    
          nFPIcnt=nFPIcnt+1
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
        if(nFPIstr .gt. 0) write(lunout,371) nFPIstr
371     format('#     Total of ',i5,' FPI strings')
        if(nFPIcir .gt. 0) write(lunout,372) nFPIcir
372     format('#     Total of ',i5,' FPI rings, free')
        if(nFPIcnt .gt. 0) write(lunout,373) nFPIcnt
373     format('#     Total of ',i5,' FPI rings with fixed centroid')

c       ---------------------------------------------------------------                   
c       ---------------------------------------------------------------                   


      rcutmax=rcut(1)
      if(rcut(2) .gt. rcutmax) rcutmax=rcut(2)
      if(rcut(3) .gt. rcutmax) rcutmax=rcut(3)
c                                                                                         
c                                                                                         
cPres                                                                                     
c  Scale initial coordinates and velocities:                                                             
      axinv=1./ax
      do i=1,3*natoms
         ra(i)=ra(i)*axinv
         va(i)=va(i)*axinv
      enddo

c  Scale other input quantities with units of length or velocity:
      if(Ideposit) then
         dz0vd=dz0vd*axinv
         v0vd=v0vd*axinv
      endif

cc      write(lunout,*) '   SETUP:    After scaling:'                                     
cc      do i=1,natoms                                                                     
cc       write(lunout,900)  i,va(3*(i-1)+1),va(3*(i-1)+2),va(3*(i-1)+3)                   
cc      enddo                                                                             
cc900   format('   atom = ',i5,'  vx,vy,vz = ',3g12.4)                                    

c   Zero piston velocity if this is a constant volume calculation:                        
       if(.not. ipres) volv=0

c   Calculate mesh size for temperature grid:                                             
       if(ITEMPGR) then
         tgrlx=ax/float(ngrlx)
         tgrly=ay/float(ngrly)
         tgrlz=az/float(ngrlz)
c                     tgrl is the mesh size of temperature grid.                          
       endif

c  Internally the temperature is stored as 0.5*D*T for simplicity (D is dimension)        
      ftemp=2./3.
      if(az .le. 0.99e-06) ftemp = 1.0
cERB      textl=textl/ftemp
cERB      temmas=temmas/ftemp
c                                                                                         
c   COMPUTE STPSZ/(2.*M) AND ITS INVERSE:                                                 
      DO 10 I=1,NATYPE
        sscrh(I)=STPSZ/(2.*AMASS(I))
10      continue
c                                                                                         
c   INITIAL MASSIVE COLLISION, IF REQUESTED:                                              
cERB      IF(INCOL) call masscol(temmas,pkinet,ekindiff)
c                                                                                         
c                   --------------------------------                                      
c                                                                                         
c   Check whether the simulation box is large enough:                                     
       if(axhalf .lt. rskinmax .or. ayhalf .lt. rskinmax) go to 998
       if(az .gt. 0.99e-06 .and. azhalf .lt. rskinmax) go to 997
cskewed       if(a2sina .lt. 2*rc .or. a1*sinalp .lt. 2*rc) go to 998                     
c                                                                                         
c-------------------------------------------------------------------------                                      
c   open more output files:                                                               
      if(IDUMP) open(luncm,file='cm.con')
      if(INSTA) open(lunins,file='ins.dat')
      if(ITEMPGR) open(luntgr,file='tgr.dat')
      if(IAVECOO) open(lunaco,file='cave.con')
      if(IVISC) open(lunvis,file='viscos.dat')

c Modified by Fer.
      if ( Iattraj ) then
        open(lunattr,file='outattraj.dat')
c       open(87,file='outattraj.xyz')
        open(89,FILE='outattraj.bin',FORM='UNFORMATTED')
c Added by Fer
        write(89) natms(2)
        write(89) ax, ay, az
      end if

      if(Iqkmin) open(lunqckmin,file='outqckmin.dat')
cc      if(nFPI .gt. 0)  open(luncoFPI,file='coFPI.con')
      if(nFPI .gt. 0)  open(lunoutFPI,file='outFPI.dat')

c              -----------------------------------------
                                                                                         
c    Write header for quickmin file (quickmin.dat):               
      if(Iqkmin) then
        WRITE(lunqckmin,240) nstatwr
240     FORMAT('# '/'#'/'#    MOLECULAR DYNAMICS   3-D EAMe  '/
     +           '#       Info on Quick Minimization:'/
     +           '#       output every  ',i5,' time steps.')

        write(lunqckmin,246)
246     format('# '/
     +     '#   Time:       RMS F:      Max|F|:  MaxFatom: ',
     +     '   x:       y:      z:        vx:     ',
     +     ' vy:      vz:        fx:        fy:      fz: ',
     +     /'#')

cc            ind=3*(indatwr-1)
cc            write(lunqckmin,241) 
cc     +                sqrt(fmag2/Natoms)*ax,fmax*ax,iatmax+1,
cc     +                Ra(ind+1)*ax,Ra(ind+2)*ax,Ra(ind+3)*ax,
cc     +                Va(ind+1)*ax,Va(ind+2)*ax,Va(ind+3)*ax,
cc     +                Fa(ind+1)*ax,Fa(ind+2)*ax,Fa(ind+3)*ax
241     format('   Initially:',g11.4,g12.4,i5,'   ',3f9.4,' ',3f9.4,' ',
     +          3f11.7)
      endif
c           close `if(Iqkmin) ... '


c  ------------------------------------------------------------------------------     
c  ------------------------------------------------------------------------------     

c   Set things up for the first time step:
      indupd=1
c          make sure a list of interactions is created now.                               
      
c      print *,'Calling FORCE'      
                                                                                   
      CALL FORCE(eQM,dEQM,FAOUT,ETOUT,natm,QPOLEOUT)
      
CAA Write out the inital forces
      open(421,file='init_forces.dat',status='replace',action='write')
      write(421,'(3g22.12)') (fa(3*i-2:3*i),i=1,natoms)
      close(421)
CAA

      if (1 .lt. 0) then

         open(87, file='coord.i',status='unknown')
         open(97, file='force.i',status='unknown')

         write(87, '(A)') '=========================================='
         do i = 1, 288
            write(87, '(3f16.9,2i6)') ra(3*(i-1)+1)*ax, ra(3*(i-1)+2)*ax
     $           ,ra(3*(i-1)+3)*ax, 0, i
         end do
         write(87, '(A)') '=========================================='
         do i = 1, 288
            write(87, '(3f16.9,2i6)') va(3*(i-1)+1)*ax, va(3*(i-1)+2)*ax
     $           ,va(3*(i-1)+3)*ax, 0, i
         end do
         
         write(97, '(A)') '=========================================='
         do i = 1, 288
            write(97, '(3f16.9,2i6)') fa(3*(i-1)+1), fa(3*(i-1)+2), fa(3
     $           *(i-1)+3)
         end do
         close(87)
         close(97)

         stop
      end if


c     Zero the force on rigid (permafrost) atoms:                                         
      call rigidatoms(dummy)
c              sero the velocity of rigid atoms.                                          
c                                                                                         
      if(naperm .gt. 0) WRITE(lunout,313) naperm
313   format('#'/
     +       '#   naperm:  Total number of rigid (permafrost) atoms = ',
     +          i10/'#')
      if(naconstr .gt. 0) WRITE(lunout,316) naconstr
316   format('#'/
     +  '#  naconstr:  Total number of constrained atoms = ',i10/'#')
c                                                                                         
c     MD2D operation:                                                                     
          if(az .le. 0.99e-06) call md2dop(dummy)
c                zero the velocity in z direction if this is a 2D simulation.             

cc      write(lunout,*) '   SETUP:    After calling force:'                               
cc      do i=1,natoms                                                                     
cc       write(lunout,901)  i,fa(3*(i-1)+1)*ax,fa(3*(i-1)+2)*ax,fa(3*(i-1)+3)*ax          
cc901   format('   atom = ',i5,'  fx,fy,fz = ',3g12.4)                                    
cc      enddo                                                                             

      Jstatwr=Nstatwr-1
cERB      if (Iqkmin)  then
cERB           write(lunqckmin,*) '   From Setup:'
cERB           call QUICKMIN
cERB      endif
c               no need to fix the velocity at this point.
c                   --------------------------------                                      

      CALL KINET
c                                                                                         
      REFKIN=TTOT
      pkinet=ttot

      if(natoms .eq. naperm) then
         write(lunout,*)' ERROR in Setup: all atoms are permafrost'
         write(lunout,*)'   natoms .eq. naperm  = ',natoms, naperm   
         stop
      endif

      temp0=ftemp*refkin/float(natoms-naperm)
      REFPOT=UTOT
      REFENG=REFKIN+REFPOT
c     ASMUS GET UTOT  
      EPOTOUT=UTOT


cPres                                                                                     
      refpre=(pkinet-0.5*virTOT)*ftemp/volm
      refpis=0.5*pmass*volv**2
cc     +                                 +pextl*volm                                      
c                  this is the total energy of the piston initially.                      

cd      write(lunout,340) time,pkinet,virTOT,volm
cd340   format('  at time= ',f8.3,'  pkinet= ',g14.6,'  virTOT= ',g14.6,
cd     +       '  volm = ',g12.4)

      refvol=volm
      refent=refeng+pextl*volm
      if(ipres .and. pmass .lt. 0.1e-30) then
         write(lunout,*) ' ERROR SETUP:  pmass must be nonzero, but is',
     +pmass
         stop
      endif
cc      if(ipres) volapp=stpsz*(refpre-pextl)/pmass                                       
c            VOLAPP is the approximate time derivative of the volume at the               
c            next timestep, t=time0+stpsz                                                 
c                                                                                         
c                 ---------------------------------------                                 
c                                                                                         
c-------------------------------------------------------------------------                                      

c  Open output files and write headers and initial values:                                

c    Write header for averages:                                                           
      write(lunout,350)
350   format('#'/'#  Averages: ')
c                                                                                         
c              -----------------------------------------
                                                                                         
c    Write header for INStantaneous values:                                               
      if(INSTA) then
        WRITE(lunins,192) NINST
192     FORMAT('# '/'#'/'#  MOLECULAR DYNAMICS   3-D   '/
     +           '#    Instantaneous values every ',i4,' steps:'/'#')
cd        write(lunins,11) (header(ii),ii=1,80)                                           
cd11      format('#',78a1)                                                                
        if(.not. IPRES) write(lunins,351) time0,temp0,refpre,
     +                                   refpot,refeng,refeng
351     format(/'#'/'#',
     +     '        Time:          Temprtr:         Pressure:',
     +     '        Potntl.En:',
     +     '       Tot.Energy:         Econst:'/
     +          '#',g22.14,5g22.14)
cc        write(lunins,*) '   refpis = ',refpis                                           
        if(IPRES) write(lunins,352) time0,temp0,refpre,refvol,
     +                              refpot,refeng,refeng+refpis
352     format('# '/'#'/
     +     '#       Time:              Temprtr:            Pressure:',
     +     '               Volume:            Potntl.En:',
     +     '           Tot.Energy:             Econst:'/
     +          '#',7g22.14)
      endif

c              -----------------------------------------
                                                                                         
c    Write header for viscosity and thermal conductivity data (viscos.dat):               
      if(IVISC) then
        WRITE(lunvis,193) NVISC,refvol
193     FORMAT('# '/'#'/'# MOLECULAR DYNAMICS   3-D   '/
     +           '#    Instantaneous values every ',i4,' steps:'/
     +           '#    Volume = ',g25.15)

cd        write(lunvis,11) (header(ii),ii=1,80)                                           
cd11      format('#',78a1)                                                                
        if(IPRES .or. IRATE) write(lunvis,*)
     +         ' WARNING: this is not a microcanonical run'
        write(lunvis,357)
357     format('# '/'#'/
     +     '#  Five lines of data are printed at each data point:'/
     +     '   # of pt:    Time:        Temprtr:     Pressure:    ',
     +     ' Energyperatom:'/
     +     '#         1/V SUM{xPx}:    1/V SUM{yPy}:  1/V SUM{zPz}:'/
     +     '#         1/V SUM{xde}:    1/V SUM{yde}:  1/V SUM{zde}:'/
     +     '#         1/V SUM{xPy}:    1/V SUM{yPx}:  1/V SUM{xPz}:'/
     +     '#         1/V SUM{zPx}:    1/V SUM{yPz}:  1/V SUM{zPy}:')

      endif
c           close `if(IVISC) ... '

c              -----------------------------------------
                                                                     
c    Write header for TGR (temperature grid) values:                                      
      if(ITEMPGR) write(luntgr,331) tgrlx,tgrly,tgrlz,ngrlx,ngrly,ngrlz,
     +        maxgrx,maxgry,maxgrz
331   format('#   Temperature grid:       Centered at (0,0,0)  '/
     +       '#    bin  10,10,10,  has the point (0,0,0) in its center'/
     +       '#   Mesh size here is: '/
     +       '#        tgrlx = ',f7.3,'  tgrly =',f7.3,'  tgrlz =',f7.3/
     +       '#                      ngrlx,ngrly,ngrlz = ',3i5/
     +       '#                     maxgrx,maxgry,maxgrz = ',3i5)
c                                                                                         
c              -----------------------------------------
                                                                                         
c    Write header for trajectory and force data (attraj.dat):               
      if(Iattraj) then
        WRITE(lunattr,230) indatwr,nstatwr
230     FORMAT('# '/'#'/'#    MOLECULAR DYNAMICS   3-D EAMe  '/
     +           '#       Trajectory of atom number:',i5/
     +           '#       output every  ',i5,' steps.')

cd        write(lunattr,11) (header(ii),ii=1,80)                                           
cd11      format('#',78a1)                                                                
        write(lunattr,231)
231     format('# '/
     +     '#   Time:         x:         y:       z:        vx:     ',
     +     ' vy:      vz:        fx:        fy:      fz:      Poten:',
     +     /'#')
         ind=3*(indatwr-1)
         write(lunattr,232) time0,
     +         ra(ind+1)*ax,ra(ind+2)*ax,ra(ind+3)*ax,
     +         va(ind+1)*ax,va(ind+2)*ax,va(ind+3)*ax,
     +         fa(ind+1)*ax,fa(ind+2)*ax,fa(ind+3)*ax,
cc     +         potenpat(indatwr),
     +         utot
232      format(' ',f10.5,'  ',3f9.4,'  ',3f9.4,'  ',3f10.5,2g22.12)
      endif
c           close `if(Iattraj) ... '

c              -----------------------------------------
                                                                                         
c    Write header for FPI data (outFPI.dat):               
      if(nFPI .gt. 0) then
        WRITE(lunoutFPI,236) nFPI,NGRP
236     FORMAT('#'/'#'/'#    MOLECULAR DYNAMICS  3-D EAMe Info on FPI'/
     +           '#       Total number of FPI chains:',i5/
     +           '#       output every  ',i5,' steps.')

        write(lunoutFPI,237)
237     format('# '/
     +     '#   Time:       avex:        avey:      avez:       ',  
     +     ' stddev:       avefx:      avefy:     avefz:        Etot:')

      endif
c           close `if(nFPI .gt. 0) ... '

                                                                                        
c  ------------------------------------------------------------------------------     
c  Write header for averages:                                                             
      pistenergy=0.5*pmass*volv**2
     +                                 +pextl*volm
c     if(.not. IPRES) write(lunout,380) time0,temp0,refpre,refpot,
c    +                                  refeng,refeng
c Modified by Fer:
c Better looking output
      if (.not. IPRES) then
        write(lunout,FMT='(A,A)')
     &'         Time     Temp        Pres       Pot. En.',
     &'      Tot. En.      Const. En.'
        write(lunout,FMT='(A,F7.2,2E12.5,3F14.8)')
     &        ' Stp: ',
     &        time0,
     &        temp0, refpre,
     &        refpot, refeng, refeng
      end if

      if(IPRES) write(lunout,381) time0,temp0,refpre,refvol,refpot,
     +                                  refeng,refeng+pistenergy
380   format('#'/'#'/'#',
     +         '    Time:   Temprtr:      Pressure:     Potntl.En:',
     +         '    Tot.Energy:       Econst: '/'#',
     +         f8.3,6g15.7)
381   format('#'/'#'/'#',
     +       '    Time:   Temprtr:      Pressure:        Volume: ',
     +       '    Potntl.En:',
     +       '    Tot.Energy:      Econst: '/'#',
     +       f8.3,6g15.7)
c                                                                                         
c                                                                                         
      RETURN

c-------------------------------------------------------------------------------
c Error messages:

998   continue
      write(lunout,9998) axhalf,ayhalf,rskinmax
9998  format(/'   ERROR in SETUP 998:  ax or ay TOO SMALL, ',
     +        ' axhalf,ayhalf,rskinmax =',
     +         4g12.4)
      stop
997   continue
      write(lunout,9997) az,azhalf,rskinmax
9997  format(/'   ERROR in SETUP 997:  az TOO SMALL,   az = ',g12.4,
     +          '  azhalf=',g12.4,'  rskinmax =',g12.4)
      stop
996   continue
      write(lunout,9996) ip,iperm(ip),va(ind+1),va(ind+2),va(ind+3)
9996  format(/'  ERROR 996 in SETUP: permafrost atom # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      stop
995   continue
      write(lunout,9995) ip,iFPI(ip),va(ind+1),va(ind+2),va(ind+3)
9995  format(/'  ERROR 995 in SETUP: last atom in FPI # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      stop
c994   continue
c      write(lunout,9994) nsumFPI,itagFPI(ich),itag(ka+indst)
c9994  format(/'  ERROR 994 in SETUP:  FPI images not adjacent ',
c     +        '  nsumFPI = ',i4,' itagFPI(ich) = ',i4,
c     +        '  itag(ka+indst) = ',i4)
c      stop
993   continue
      write(lunout,9993) iconstr(ip),va(ind+1),va(ind+2),va(ind+3)
9993  format(/'  ERROR 993 in SETUP:  constrained atom ',i5/
     +        '           has illegal velocity,  vx,vy,vz = ',3g15.5)
      stop
990   continue
      write(lunout,9990) i,itag(i)
9990  format(/'  ERROR 990 in SETUP: illegal tag on atom ',i4,':'/
cc     +        '     RA = ',3g15.6,
     +        '    itag = ',i10)
989   continue
      write(lunout,9989) ip,iFPI(ip),va(ind+1),va(ind+2),va(ind+3)
9989  format(/'  ERROR 989 in SETUP: first atom in FPI # ',i4/
     +        '     which is atom # ',i4,' has nonzero velocity:'/
     +        '     vx,vy,vz = ',3g12.4)
      stop
      END
 

