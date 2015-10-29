c   3dme:                                                                                 
c   3dmd:                                                                                 
C     THIS ROUTINE PERFORMS THE DYNAMICS FOR ntotstps TIMESTEPS                           
c                                                                                         
      SUBROUTINE LOOP
c                                                                                         
      implicit real*8 (a-h,o-z)
c                                                                                         
c                                                                                         
      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comoutcntrl.cmn'
      include '../commonblks/comtgr.cmn'
      include '../commonblks/compotent.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comenperat.cmn'
      include '../commonblks/comdeposit.cmn'
      include '../commonblks/comintlis.cmn'
      include '../commonblks/constraints.cmn'


      COMMON /ISEEDS/ ISEED1
      COMMON /REFENC/ REFPRE,REFVOL,REFPOT,REFKIN,
     +                REFENT,REFENG,temp0,refpis
      common /cputim/cput0,cputas,cputil,cputig,cpudel,cputr2,cputbl,
     +               cpudef,cputpb,cpupbs,cpucpr,cpu372,cpu370,cpu373
     +              ,cpu380,cpu381,cpu382,cpu383,cpu375,cpu377,cpu376
     +              ,cpu3761
      real*4 cputf2(2), cputbl2(2)

      dimension sumtgr(MAXTGR,MAXTGR,MAXTGR),
     +     avengr(MAXTGR,MAXTGR,MAXTGR),rdispl(MAXCOO),
     +     r2disp(MAXATOMS),natmsold(MAXCOMP),itagst(MAXATOMS),
     +     aveposFPI(MAXFPICOO),avedistFPI(MAXFPICOO),
     +     aveforcFPI(MAXFPICOO),FPIprevim(MAXFPICOO),displim(MAXFPI),
     +     relenimage(MAXFPI),ioldFPIat(MAXATOMS)

cERB Temp. grid ---(start)---
      integer ngr(MAXTGR,MAXTGR,MAXTGR)
cERB Temp. grid ---(end)---
      integer iatfla(MAXATOMS)
c            iatfla -     list of atoms that fly away.                                    
      logical iwait,IWVISC

cG77      type TB_TYPE                      !
cG77        SEQUENCE                        !All this has to be included
cG77        REAL(4) USRTIME                 !at the beginning of the
cG77        REAL(4) SYSTIME                 !code to define the structure
cG77      END TYPE                          !that the function etime_
cG77                                        !likes.
cG77      TYPE(TB_TYPE) ETIME_SRTUCT        !

c--------------------------------------------------------------------H2O
c               Lines added for simulation of H2O

      real*8 r1(3), r2(3), r3(3), v1(3), v2(3), v3(3), box(3)
      real*8 RAold(maxCoo), VAold(maxCoo), FAold(maxCoo)
      real*8 InertiaMoment(3), angularVel(3), vCM(3), pi, mMolecule
      real*8 zCoord(maxCoo), zCoord2(maxCoo)
      integer chainLength, histogram(1000), binNumber
      integer iVelAtms(maxCoo/3), nWaitVel
      logical H2O, iSaveVel
      character*80 dummyCH

C Added by Fer.
      integer*4 Pt_i, Pt_i1, Pt_i2
      real*8    Ri(3), Ri1(3), Ri2(3)
      integer*4 i, p

c     recognize if we are dealing with water

      pi = 3.14159265358979312
      if((NATYPE .eq. 2) .and. (NATMS(1) .eq. 2*NATMS(2))) then
         H2O = .TRUE.
c        compute the inertia tensor
         mMolecule = 2.d0 * amass(1) + amass(2)
         Ro = potpar(14, 1)
         theta0 = potpar(15, 1)
         h = Ro * cos(theta0/2.d0 * pi/180.d0)
         a = 2.d0 * Ro * sin(theta0/2.d0 * pi/180.d0)

         InertiaMoment(1) = 2.d0 * amass(1) * amass(2) * h*h / mMolecule
         InertiaMoment(2) = amass(1) * a*a / 2.d0
         InertiaMoment(3) = InertiaMoment(1) + InertiaMoment(2)
         nOxygens = nAtms(2)
         nHydrogens = nAtms(1)
         nMolecules = nOxygens

         box(1) = ax
         box(2) = ay
         box(3) = az
         if (IrigidMolecules) then 
            print '(A)', 'Simulating RIGID water molecules'
         else
            print '(A)', 'Simulating FLEXIBLE water molecules'
         end if
      else
         H2O = .FALSE.
      end if
c--------------------------------------------------------------------H2O

c The initial massive collisions is done in this file now. ERB

      IF(INCOL) then
         print '(A)', ' MAKE INITIAL MASSIVE COLLISIONS'
c --------
c For H2O:
         IF (H2O .and. IrigidMolecules) then
            call masscolH2O(temmas, pkinet, ekindiff, InertiaMoment, 
     $           mMolecule)
            if (genConstraints) call projectForce
c ---------
         else
            call masscol(temmas, pkinet, ekindiff)
         end if
      end if

c-----------------------------------------------------------------------
       write(6,*) '    Running the loop over time steps ... '

       ftemp = 2./3.
       if(az .le. 0.99e-06) ftemp= 1.0

c     H2O: 6 degrees of freedom for every 6 atoms (2 hydrogens and 1 oxygen)
       if (IrigidMolecules) ftemp = 1.0

c-----------------------------------------------------------------------

cG77      call second(cputbl)
c666      call etime(cputbl2)
c
          volm0=ax0*ay0*az0

          pressure=refpre

          axhalf=0.5*ax
          ayhalf=0.5*ay
          azhalf=0.5*az

          axinv=1./ax
          axsc=1.0
          aysc=ay/ax
          azsc=az/ax

          axhalfsc=0.5
          ayhalfsc=0.5*aysc
          azhalfsc=0.5*azsc


C   INITIALIZE ACCUMULATORS FOR QUANTITIES TO BE AVERAGED                               
C     
      maxgrx=(ngrlx-1)/2
      maxgry=(ngrly-1)/2
      maxgrz=(ngrlz-1)/2
      do 90 ix=1,MAXTGR
        do 91 iy=1,MAXTGR
          do 92 iz=1,MAXTGR
             sumtgr(ix,iy,iz)=0.
             avengr(ix,iy,iz)=0.
 92       continue
 91    continue
 90   continue
      

      do 120 i=1,MAXCOO
        sumRA(i)=0.0
        rdispl(i)=0.0
        raNoPBC(i)=ra(i)
120     continue
      nstupd=0

      do if=1,MAXFPI
         ind=3*(if-1)
         aveposFPI(ind+1)=0.0
         aveposFPI(ind+2)=0.0
         aveposFPI(ind+3)=0.0
         avedistFPI(ind+1)=0.0
         avedistFPI(ind+2)=0.0
         avedistFPI(ind+3)=0.0
         aveforcFPI(ind+1)=0.0
         aveforcFPI(ind+2)=0.0
         aveforcFPI(ind+3)=0.0
      enddo

      SUMTEM=0.0
      SSQTEM=0.0
      SUMPOT=0.0
      SSQPOT=0.0
      SUMENG=0.0
      SSQENG=0.0
      SUMCNS=0.0
      SSQCNS=0.0
      sumpre=0.0
      ssqpre=0.0
      sumvol=0.0
      ssqvol=0.0
c                                                                                         
      v2perm=0.0
c                                                                                         
      nsal=0
      saltem=0.0
      salpot=0.0
      saleng=0.0
      salcns=0.0
      salpre=0.0
      salvol=0.0
      RANKIN=0.0
      ranpis=0.0
      disp2addper=0
      dispmaxaddper=0
      naddper=0
      distspat=0
      sdistspat=0
      s2distspat=0

c  Surface, VD features:
      vdkin=0.0
      remkin=0.0
      flwkin=0.0
      nflawy=0
c        count the number of atoms that fly away from the surface.

C                                                                                         
C   INITIALIZE COUNTERS:                                                                
C         
      JGRP=0
      JUSR=0
      JINST=0
      JMASV=0
      JDUMP=0
      JVISC=0
      Jstatwr=0
      NVISCPT=0
      javecoo=0
      navecoo=0
      nstochcoll=0
      nmasscoll=0
      IWAIT=.true.
      IWVISC=.true.
      TIME=time0
      jvdep=ntotvd-1
c        set the counter for Atom Deposition so that an atom will be deposited            
c        right away.                                                                      
C                                                                                         
      icntud=0
      indaddper=0

c-----------------------------------------------------------------------
      if (.false.) then
         open(87, file='coord.i',status='unknown')
         open(97, file='force.i',status='unknown')
      end if
c-----------------------------------------------------------------------
c                        MD loop begins
c-----------------------------------------------------------------------
c  Begin the big loop over timesteps:                                                     
cTimming
cTimming      time0 = etime_(ETIME_STRUCT)
cTimming      time1 = 0.
cTimming      time2 = 0.
cTimming      time3 = 0.
cTimming      time4 = 0.
cTimming      time5 = 0.
cTimming      time6 = 0.
cTimming      time7 = 0.
cTimming

      iterations = -1

cHisto      z0 = -10.
cHisto      dz = 0.02

c------ Saving velocities ----
      iSaveVel = .FALSE.
c      iSaveVel = .TRUE.
      if (iSaveVel) then
         open (89, file = '/scratch/erb/cmv.con')
c         open (89, file = 'cmv.con')
         open (90, file = 'velAtms.dat')
         read (90,*) dummyCH
         do i = 1, natoms
            read (90,*, err=11, end=11) iVelAtms(i)
         end do
 11      close(90)
         velsave = -1
         nWaitVel = 4
         nMolSave = i-1
c$$$         print *, nMolSave
c$$$         do i = 1, nMolSave
c$$$            print *, iVelAtms(i)
c$$$         end do
      end if
c-----------------------------

c-----Generalized Constraints-------
c$$$      if (genConstraints) then
c$$$         vDotA = 0.d0
c$$$         do i = 1, 3*nAtoms
c$$$            vDotA = vDotA + VA(i)*constraint(i)
c$$$         end do
c$$$         do i = 1, 3*nAtoms
c$$$            va(i) = va(i) - vDotA * constraint(i)
c$$$         end do
c$$$      end if

      if (genConstraints) call projectForce
         
      If ( Iattraj ) then

C Added by Fer.
C Dump important stuff
        write(89) time, uTot, (ax*Ra(i),ax*Va(i),i=1,9*nMolecules)

      end if

c*********** MD loop ************
      DO 100 itime=1,ntotstps
c---------------------------------
c---- Dump position of molecule 1 -----
c         write(88,'(9f15.8)') ra(1), ra(2), ra(3), ra(4), ra(5), ra(6),
c     $        ra(3*natms(1)+1), ra(3*natms(1)+2), ra(3*natms(1)+3)
c--------------------------------------



c------ Saving velocities ----
         if (iSaveVel) then
            velsave = velsave + 1
            if (velsave .eq. nWaitVel) then
               call velout(89, ax, iVelAtms, nMolSave)
               velsave = 0
            end if
         end if
c-----------------------------
         

cTimming
cTimming      time1 = etime_(ETIME_STRUCT)
cTimming      write (*, '(A,f15.10)') 'Total run time: ', time1-time0
cTimming      write (*, '(A,f15.10)') 'Miscellaneous times: ', 
cTimming     $     (time1-time0) - (time3-time2) - (time5-time4) - (time7-time6)
cTimming      write (*,*)
cTimming      write (*,*)
cTimming

c     Update the HISTOGRAM
c$$$         do j = nHydrogens + 1, nHydrogens+nOxygens
c$$$            z = ra(3*j) * ax
c$$$            binNumber = int((z-z0)/dz) + 1
c$$$            histogram(binNumber) = histogram(binNumber) + 1
c$$$         end do

c$$$         iterations = iterations + 1
c$$$         timeFinal = etime_(ETIME_STRUCT)
c$$$
c$$$         print '(f15.10)', iterations, timeFinal-timeInitial
c$$$         timeInitial = timeFinal

c Save z displacement of movile oxygens         

c$$$         do iatom = nHydrogens+1, nHydrogens+nOxygens
c$$$            zCoord2(iatom-nHydrogens) = zCoord2(iatom-nHydrogens) + 
c$$$     $           ra(3*iatom) * ra(3*iatom)
c$$$            zCoord(iatom-nHydrogens) = zCoord(iatom-nHydrogens) + 
c$$$     $           ra(3*iatom)
c$$$         end do

c     Save previous coordinates, velocities and forces. If doing
c     quickmin and overshoots, the program will come back to the
c     previous configuration. (BUT IT IS NOT WORKING YET => comments)
         do i = 1, 3*natoms
c            FAold(i) = FA(i)
            RAold(i) = RA(i)
c            VAold(i) = VA(i)
         end do
                                                                                         
C      UPDATE THE COORDINATES: (NB: the coordinates are all scaled)                       
       ISTR=1
       DO 105 ITYPE=1,NATYPE
        DO 104 IATOM=ISTR,ISTR+3*NATMS(ITYPE)-1
104        FA(IATOM)=SSCRH(ITYPE)*FA(IATOM)
c                   'sscrh' (on commonblk 'time') contains the stepsize                   
c                                   over twice the mass.                                  
105     ISTR=ISTR+3*NATMS(ITYPE)

cc       go to 9000                                                                       

cPres:                                                                                    
      fact=1.0
      if(IPRES)  fact=1.0-((stpsz*volv)/(3.*volm))

        do I=1,3*NATOMS
           VA(I)=fact*VA(I)+FA(I)
           RA(I)=RA(I)+(STPSZ*VA(I))
c                   the scaled coordinates have now been updated to new time.             
           raNoPBC(i)=raNoPBC(i)+(STPSZ*VA(I))
           sumRA(i)=sumRA(i)+RA(i)
c               calculate the average coordinates of the atoms.                           
           rdispl(i)=rdispl(i)+(STPSZ*VA(I))
c               calculate the total displacement since last update of list.               
        enddo
         navecoo=navecoo+1
c               count the number of steps included in sumRA.                              

cc        write(lunout,*)'#    LOOP:   coordinates updated.'                                  
cc      do i=1,natoms                                                                     
cc         i=1
cc         write(lunout,901)  i,ra(3*(i-1)+1)*ax,ra(3*(i-1)+2)*ax,
cc     +                      ra(3*(i-1)+3)*ax          
901      format('   atom = ',i5,'  RAx,RAy,RAz = ',3g12.4)                                    
cc      enddo                                                                             
cTimming
cTimming      time2 = etime_(ETIME_STRUCT)
cTimming
c
c-----------------------------------------------------------------------
c     Use the RATTLE algorithm to simulate rigid molecules.
c     Here we call the first part, i.e. the corrections to r(t+dt) and
c     corrections for v(t+dt/2)
c         print '(A, $)', 'Are molecules rigid? '
c         print *, IrigidMolecules

         if (.false.) then
         print *, 'to shaker1'

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
c$$$         close(87)
c$$$         close(97)
         end if
         
C Debug
C     write(10,FMT='(1X,A)') 'Before shaker1'
C     write(10,FMT='(1X,A)') 'RA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (RA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'VA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (VA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'FA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (FA(3*(i-1)+p),p=1,3)
C     end do

         if (H2O .and. IrigidMolecules) call shaker1(RAold)

C Debug
C     write(10,FMT='(1X,A)') 'After shaker1'
C     write(10,FMT='(1X,A)') 'RA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (RA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'VA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (VA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'FA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (FA(3*(i-1)+p),p=1,3)
C     end do

c-----------------------------------------------------------------------
cTimming
cTimming      time3 = etime_(ETIME_STRUCT)
cTimming      write (*, '(A,f15.10)') 'Time spent in Shaker1: ', time3-time2
cTimming
c                                                                                         
cPres:                                                                                    
      if(IPRES) then
c      update the volume:                                                                 

c       Estimate the time derivative of volume for the current step using                 
c       values of volv and pressure from the last step:                                   
          volapp=volv+stpsz*(pressure-pextl)/pmass

cc          write(lunout,342) volapp                                                      
cc342       format('                volapp =  ',g12.4)                                    

          volm=volm+0.5*stpsz*(volv+volapp)
c           update scalefactor for distances:  (use ax as scale factor)                   
          if(volm .lt. 0.) then
           write(lunout,*) '**  ERROR in LOOP:   Volm < 0, volm = ',volm
           stop
        endif
          scale=(volm/volm0)**(1./3.)
          ax=ax0*scale
          ay=ay0*scale
          az=az0*scale
c2d	  if(ind2D) az=1.0                                                                    
          axhalf=0.5*ax
          ayhalf=0.5*ay
          azhalf=0.5*az

          axinv=1./ax
          axsc=1.0
          aysc=ay*axinv
          azsc=az*axinv
          axhalfsc=0.5
          ayhalfsc=0.5*aysc
          azhalfsc=0.5*azsc


cc          write(lunout,*)'#   LOOP:    volume updated, volm=',volm                       

      endif
c                                                                                         


c               ----------------------------------------                                  
c                                                                                         
c  PBC, apply periodic boundary conditions:                                               
      if(abs(cosalp) .lt. 0.1e-10) then
c                                                                                         
      aynavecoo=aysc*navecoo
      axnavecoo=axsc*navecoo
      aznavecoo=azsc*navecoo
c    Rectangular simulation cell:                                                         
      DO 110 I=1,NATOMS
         ind=3*(i-1)
c      first the y coordinate:                                                            
         IF(RA(ind+2) .GT. ayhalfsc) then
            RA(ind+2)=RA(ind+2)-aysc
            sumra(ind+2)=sumRA(ind+2)-aynavecoo
         endif
         IF(RA(ind+2) .LE. -ayhalfsc) then
            RA(ind+2)=RA(ind+2)+aysc
            sumRA(ind+2)=sumRA(ind+2)+aynavecoo
         endif
c      then the x coordinate:                                                             
         IF(RA(ind+1) .GT. axhalfsc) then
            RA(ind+1)=RA(ind+1)-axsc
            sumRA(ind+1)=sumRA(ind+1)-axnavecoo
         endif
         IF(RA(ind+1) .LE. -axhalfsc) then
            RA(ind+1)=RA(ind+1)+axsc
            sumRA(ind+1)=sumRA(ind+1)+axnavecoo
         endif
c      finally the z coordinate:                                                          
         IF(RA(ind+3) .GT. azhalfsc) then
            RA(ind+3)=RA(ind+3)-azsc
            sumRA(ind+3)=sumRA(ind+3)-aznavecoo
         endif
         IF(RA(ind+3) .LE. -azhalfsc) then
            RA(ind+3)=RA(ind+3)+azsc
            sumRA(ind+3)=sumRA(ind+3)+aznavecoo
         endif
110      CONTINUE
c                                                                                         
       else
c                                                                                         
c   DEVELPMNT:     This has to be extended to 3-D                                         
c     Skewed simulation cell:                                                             
      write(lunout,*) 
     +     '#** WARNING from LOOP:  Skewed simulation cell not ready'
      write(lunout,*) '#   cosalp = ',cosalp
      DO 111 I=1,NATOMS
         ind=3*(i-1)
c      first the y coordinate:                                                            
         IF(RA(ind+2) .GT. ayhalfsc) then
               RA(ind+1)=RA(ind+1)-ay*cosalp
               RA(ind+2)=RA(ind+2)-ay*sinalp
         endif
         IF(RA(ind+2) .LE. -ayhalfsc) then
               RA(ind+1)=RA(ind+1)+ay*cosalp
               RA(ind+2)=RA(ind+2)+ay*sinalp
         endif
c      then the x coordinate:                                                             
         IF(RA(ind+1) .GT. axhalf+RA(ind+2)*cosalp/sinalp)
     +                            RA(ind+1)=RA(ind+1)-ax
         IF(RA(ind+1) .LE.-axhalf+RA(ind+2)*cosalp/sinalp)
     +                            RA(ind+1)=RA(ind+1)+ax
111      CONTINUE
       endif
c                                                                                         

cc        write(lunout,*)'#   LOOP:    Apply PBC'                                          

c                  ---------------------------------------           end PBC.                
       if (.not. H2O) then             ! SKIP neighbor list in H2O

c   NEIGHBORLIST:                                                                         
c       Check whether the list of interacting pairs should be updated.  Find              
c       the sum of the magnitude of the two largest displacements and compare             
c       with the difference in potential cutoff and skin radius:                          
         do 121 i=1,natoms
           r2disp(i)=rdispl(3*i-2)**2+rdispl(3*i-1)**2+rdispl(3*i)**2
121        continue
         dmax=0.
         dmaxm=0.
c     find the largest displacement:                                                      
         do 122 i=1,natoms
           if(r2disp(i) .le. dmax) go to 122
              dmax=r2disp(i)
              istm=i
122        continue
c     find the second largest displacement:                                               
         r2disp(istm)=0.
         do 123 i=1,natoms
           if(r2disp(i) .le. dmaxm) go to 123
              dmaxm=r2disp(i)
              istmm=i
123        continue
         summag=sqrt(dmax)+sqrt(dmaxm)
c                                                                                         
cq      write(lunout,401) dmax,istm,dmaxm,istmm,(r2disp(i),i=1,natoms)                    
cq401   format(/'   Find largest displacement:  dmax=',g12.4,' atom ',i4,                 
cq    +        ';  dmaxm=',g12.4,' atom ',i4/                                             
cq     +        ' r2disp:',(7g11.3))                                                      
c                                                                                         
c   Find wich components the two atoms belong to: (assume only two comp.)                 
         if(istm .le. natms(1) .and. istmm .le. natms(1)) indcom=1
         if(istm .gt. natms(1) .and. istmm .gt. natms(1)) indcom=2
         if(istm .le. natms(1) .and. istmm .gt. natms(1)) indcom=3
         if(istm .gt. natms(1) .and. istmm .le. natms(1)) indcom=3
c            indcom=1 for AA,  indcom=2 for BB,  indcom=3 for AB or BA.                   
c                                                                                         
c     Have to update the list if any pair of atoms has moved far enough apart:            
        if(skndpt(indcom) .lt. summag*ax) then
           indupd=1
c                 when the value of indupd is 1, then the list of interacting             
c                 pairs will be updated next time GAGAFE is called.                       
c        reset the displacement:                                                          
           do 124 i=1,3*natoms
             rdispl(i)=0.
124          continue
         endif
c                                                                                         
         end if          ! SKIP neighbor list in H2O
c      ---------------------------------------------------------- end neighb. list.                  
C                                                                                         
C   COMPUTE NEW FORCES AND NEW PARTICLE VELOCITIES:                                       
cTimming
cTimming      time4 = etime_(ETIME_STRUCT)
cTimming

C Debug
C     write(10,FMT='(1X,A)') ' Before force '
C     write(10,FMT='(1X,A)') 'RA'
C     do i=1,natoms
C       write(10,FMT='(1X,3E23.15)') (RA(3*(i-1)+p),p=1,3)
C     end do

      CALL FORCE
c$$$      if (genConstraints .and. iqkmin) call projectForce
      if (genConstraints) call projectForce
cTimming
cTimming      time5 = etime_(ETIME_STRUCT)
cTimming      write (*, '(A,f15.10)') 'Time spent in Forces: ', time5-time4
cTimming

c    finish updating particle velocities:                                                 
      ISTR=1
      fact=1.0
      if(IPRES) fact=1./(1.+stpsz*volapp/(3.*volm))
      DO 135 ITYPE=1,NATYPE
         DO IATOM=ISTR,ISTR+3*NATMS(ITYPE)-1
             VA(IATOM)=fact*(VA(IATOM)+SSCRH(ITYPE)*FA(IATOM))
         enddo
135      ISTR=ISTR+3*NATMS(ITYPE)
c
cTimming
cTimming      time6 = etime_(ETIME_STRUCT)
cTimming
c-----------------------------------------------------------------------
c     Use the 2nd part of the RATTLE algorithm to update the velocities
c     of the rigid molecules, i.e. v(t+dt).

C Debug
C     write(10,FMT='(1X,A)') 'Before shaker2'
C     write(10,FMT='(1X,A)') 'RA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (RA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'VA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (VA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'FA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (FA(3*(i-1)+p),p=1,3)
C     end do

         if (H2O .and. IrigidMolecules) call shaker2()

C Debug
C     write(10,FMT='(1X,A)') 'After shaker2'
C     write(10,FMT='(1X,A)') 'RA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (RA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'VA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (VA(3*(i-1)+p),p=1,3)
C     end do
C     write(10,FMT='(1X,A)') 'FA'
C     do i=1,natoms
C       write(10,FMT='(1X,3F20.14)') (FA(3*(i-1)+p),p=1,3)
C     end do

c-----------------------------------------------------------------------
cTimming
cTimming      time7 = etime_(ETIME_STRUCT)
cTimming      write (*, '(A,f15.10)') 'Time spent in Shaker2: ', time7-time6
cTimming


c       Zero velocity and force of rigid atoms:
          call rigidatoms(v2perm)

c      -----------------------------------------------------------------
c
C   COMPUTE KINETIC ENERGY:
C
      CALL KINET
c
      PKINET=TTOT
      PPOTEN=UTOT

cPres:
      pressure=(pkinet-0.5*virTOT)*ftemp/volm

c$$$      print '(A,f15.8)', 'Virial(main) = ',virTot
c$$$      print *


      if(IPRES) then
         volv=0.5*(volv+volapp)+(0.5*stpsz/pmass)*(pressure-pextl)
c           Now the time derivative of the volume has been updated and that               
c           completes the update of generalized coordinates and momenta for new           
c           time.                                                                         

c        Piston collisions: (piston suffers collisions at every step if IPRES)            
c                              NB!  bath temperature is needed here even if the           
c                                   atoms suffer no collisions.                           
         pistenergy=0.5*pmass*volv**2
     +                                 +pextl*volm
         if(IPCOL .and. IPRES) then
              volvdyn=volv
              call ranvel(volvstoch,dummy,dummy,pmass,textl)
c             scale:                                                                      
              volvstoch=volvstoch/ax
              volv=sqrt(1.0-plambd**2)*volvdyn+plambd*volvstoch
              ranpis=ranpis+0.5*pmass*(volvdyn**2-volv**2)

cd              write(lunout,341) volvdyn,volvstoch,volv,plambd                           
cd341           format('#  Piston collision:   volvdyn =',g12.4,' volvstoch =',            
cd     +                g12.4,' volv =',g12.4/'#        plambd = ',g12.4)                  


         endif
c                       close 'if(IPCOL .and. IPRES)'                                     

      endif

c    -----------------------------------------------------------------                    

C  STOCHASTIC COLLISIONS ON THE ATOMS:  MASSIVE AND/OR RATEWISE                           

c    First,   MASSIVE collisions:                                                         
      IF(IMASV) then
c            This condition corresponds to periodic massive collisions.                   
        JMASV=JMASV+1
c            this variable counts the number of time steps since last MASS COLL           
        IF(JMASV .eq. NMASV) then
          JMASV=0

cd        write(6,*) '#   just before masscol:    naperm = ',naperm                        

          nmasscoll=nmasscoll+1
c --------
c For H2O:
         IF (H2O .and. IrigidMolecules) then
            call masscolH2O(textl, pkinet, ekindiff, InertiaMoment, 
     $           mMolecule)
            if (genConstraints) call projectForce
c ---------
         else
            call masscol(textl, pkinet, ekindiff)
         end if
cERB          call masscol(textl,pkinet,ekindiff)
c                       The new velocities come already scaled from masscol.              
          rankin=rankin-ekindiff
c              keep track of the change in energy due to stochastic collisions.           

cPres:   The kinetic energy is changed in 'masscol', then pressure is changed:             
          pressure=(pkinet-0.5*virTOT)*ftemp/volm
c                                     the PV energy is changed as a result.               

cd        write(lunout,340) time,pkinet,virTOT,volm                                         
cd340     format('#   After masscol at time = ',f8.3,
cd     +         ' pkinet ',g14.6,' virTOT =',           
cd     +           g14.6,'  volm = ',g12.4)                                               

c       The estimate of the time derivative of volume for the next step                     
c        needs to be revised: (just as if the calculation was being started now           
c        using the reassigned velocities as starting velocities):                         
          if(IPRES) then
cd          volv=0.5*(volv+volapp)+(0.5*stpsz/pmass)*(pressure-pextl)                      
            volapp=volv+stpsz*(pressure-pextl)/pmass
          endif

cd      write(6,*)'# LOOP:   after masscol at t=',time,' ekindiff = ',ekindiff        

        endif
      endif
c            close 'if(imasv)'.                                                           

cDeposit, VD:
       if(Ideposit .and. (jvdep .eq. naftvd)) then
c            this condition corresponds to the first massive collision 
c            after a new atom was sent towards the surface in atom deposition.
c --------
c For H2O:
         IF (H2O .and. IrigidMolecules) then
            call masscolH2O(textl, pkinet, ekindiff, InertiaMoment, 
     $           mMolecule)
            if (genConstraints) call projectForce
c ---------
         else
            call masscol(textl, pkinet, ekindiff)
         end if
cERB	  call masscol(textl,pkinet,ekindiff)  
	  rankin=rankin-ekindiff 
c      Get some info on the fate of the deposited atom:
         if(natdep .ge. 1) then
           indspat=1
           do it=1,10
             if(itypvd .eq. 2) indspat=1+natms(1)
           enddo
c         PBC:
           xdist=abs(x0vd-RA(3*(indspat-1)+1))
              if(xdist .gt. axhalfsc) xdist=xdist-axsc
           ydist=abs(y0vd-RA(3*(indspat-1)+2))
              if(ydist .gt. ayhalfsc) ydist=ydist-aysc
           distspat2=(xdist)**2+(ydist)**2             
           distspat=sqrt(distspat2)*ax
           write(lunout,510) distspat
510        format('#   distance traveled by deposited atom in (x,y)=',
     +                   g12.5)
           if(distspat .gt. 2) write(lunout,512) RA(3*(indspat-1)+1),
     +                                 RA(3*(indspat-1)+1),x0vd,y0vd
512        format('#    WOW!  This atom has moved far, now (x,y) =',
     +              2g12.4/'#          but initially (x,y) = ',2g12.4)

           sdistspat=sdistspat+distspat
           s2distspat=s2distspat+distspat2
         endif
       endif       
       if(Ideposit .and. (jvdep .eq. (naftvd+nbtwmc))) then
c            this condition corresponds to second mass. coll. after atom dep.
c --------
c For H2O:
         IF (H2O .and. IrigidMolecules) then
            call masscolH2O(textl, pkinet, ekindiff, InertiaMoment, 
     $           mMolecule)
            if (genConstraints) call projectForce
c ---------
         else
            call masscol(textl, pkinet, ekindiff)
         end if
cERB	  call masscol(textl,pkinet,ekindiff)
	  rankin=rankin-ekindiff 
        endif
c             - - - - - - -  - - - - - - - - -  - - - - - - - -                           

c   Second,   INDIVIDUAL stochastic collisions:                                           


c$$$        do iii = 1, 6
c$$$           print '(i4, 3f15.8)', iii, (va(3*(iii-1)+j), j = 1, 3)
c$$$        end do
c$$$        pause

      IF(IRATE .and. nlangevin .eq. 0) then
C        CRATE is the # of collisions per timestep:                                       
      J=INT(CRATE)

      F=CRATE-J
c        When CRATE is not an integer, the number of atoms suffering a collision          
c        is not always the same.  F is used to                                            
c         decide whether there will be J or J+1 collisions at this step:                  
      IF(RANb(ISEED1,ISEED2).LT.F) J=J+1
      nstochcoll=nstochcoll+J

c    Stochastic collision with J atoms:

       ntries=0
       DO 255 I=1,J
254     continue
        ntries=ntries+1

c----H2O /*
cERB        IATOM=FLOAT(NATOMS)*RANb(ISEED1,ISEED2)+1.
c---------
c For H2O:
        if (H2O .and. IrigidMolecules) then
c          H2O: Choose a molecule instead of an atom
           IATOM=FLOAT(nOxygens)*RANb(ISEED1,ISEED2)+1.
        else
           IATOM=FLOAT(NATOMS)*RANb(ISEED1,ISEED2)+1.
        end if
c----H2O */

c   make sure the atom is in the region where temperature is controled:
        if(Ideposit) then
          if(ntries .gt. 1000) then
            write(lunout,*)
     $            'Problem in fining atoms to hit, ntries=',ntries
            stop
          endif
          if(ra(3*(IATOM-1)+3)*ax .gt. ztcontrl) go to 254
        endif

        NAT=0
c     find what type the atom is:                                                         
c-------------H2O /*
c     For H2O:
c     It is not necesary to find the type of atom if playing with rigid
c     molecules. It is the whole molecule the one that collides. ERB.
        if ((.not. H2O) .or. (.not. IrigidMolecules)) then

        DO 251 ITYPE=1,NATYPE
         NAT=NAT+NATMS(ITYPE)
         IF(IATOM.GT.NAT)GO TO 251
          ITY=ITYPE
          GO TO 253
251      CONTINUE
       write(lunout,*) 
     +    ' ERROR in LOOP:  could not find type of atom ',iatom
       stop
c                                                                                         
253     continue
        end if
c-------------H2O */

cd        write(6,*) '#   LOOP:   just before stoch coll. on atom ',iatom,                 
cd     +         '  at time = ',time                                                      
cd        write(6,*) '#         Va = ',VA(3*IATOM-2),VA(3*IATOM-1),VA(3*IATOM)           



c-------------H2O /*
c--------------------------------------------Start of part for rigid H2O
c Coupling to a heat bath using 6 degrees of freedom, the coordinates of
c the center of mass and the three angles that define the angular
c velocity. ERB

        if (IrigidMolecules) then

           chainLength = nimFPI(1)
           if (nimFPI(1) .eq. 0) chainlength = 1
           moleNumber = int((iatom-1)/chainLength) + 1
           image      = iatom - chainLength * (moleNumber - 1)
           iatom = moleNumber

           indexO = 3 * (nAtms(1) + (iatom-1) * chainLength + 
     $          (image-1))
           indexH1 = 3 * (2 * (iatom-1) * chainLength + (image-1))
           indexH2 = indexH1 + 3 * chainLength
           do ii = 1, 3
              r1(ii) = RA(indexH1 + ii)
              r2(ii) = RA(indexH2 + ii)
              r3(ii) = RA(indexO + ii)
c$$$              r1(ii) = RA(6*(iatom-1) + ii)
c$$$              r2(ii) = RA(6*(iatom-1) + 3 + ii)
c$$$              r3(ii) = RA(3*(nHydrogens + iatom - 1) + ii)
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
     $           mMolecule)

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
c$$$            do ii = 1, 3
c$$$               index = 6*(iatom-1) + ii
c$$$               indexO = 3*(iatom-1) + 3*nHydrogens + ii
c$$$               VA(index)   = v1(ii) + vCM(ii)
c$$$               VA(index+3) = v2(ii) + vCM(ii)
c$$$               VA(indexO)  = v3(ii) + vCM(ii)
c$$$            end do

c------------------------------------------------- End of part for rigid H2O
        else
c     Non-constrained dynamics => draw velocities at random
           CALL RANVEL(vx,vy,vz,AMASS(ITY),TEXTL)
           VA(3*IATOM-2) = vx * axinv
           VA(3*IATOM-1) = vy * axinv
           VA(3*IATOM)   = vz * axinv
        end if
255   CONTINUE
c------------------H2O */


c    Rigid atoms:                                                                         
      call rigidatoms(dummy)
c              zero the velocity of rigid atoms.                                          
c2     2DMD operation:                                                                    
c2          if(IND2D) call md2dop(dummy)                                                  
c2                zero the velocity in z direction if this is a 2D simulation.            

      CALL KINET

      PKNEW=TTOT

      ekdiff=(PKINET-PKNEW)

cd      write(6,*) '#        time = ',time,'   ekdiff = ',ekdiff                           


c$$$        do iii = 1, 6
c$$$           print '(i4, 3f15.8)', iii, (va(3*(iii-1)+j), j = 1, 3)
c$$$        end do
c$$$        pause


      RANKIN=RANKIN+(PKINET-PKNEW)
c          this is the change in kinetic energy due to the stoch. coll.                   
      PKINET=PKNEW

      pressure=(pkinet-0.5*virTOT)*ftemp/volm

      if(IPRES) then
cd         volv=0.5*(volv+volapp)+(0.5*stpsz/pmass)*(pressure-pextl)                      
c        Estimate of the time derivative of volume for the next step                      
c        needs to be revised: (just as if the calculation was being started now           
c        using the reassigned velocities as starting velocities):                         
         volapp=volv+stpsz*(pressure-pextl)/pmass
      endif
C                                                                                         
      endif
c           close 'if(irate)' 
      if (genConstraints) call projectForce
c            - - - - - - -  - - - - - - - - -  - - - - - - - - 
                          

c   Third,   Langevin stochastic collisions:                                           
      IF(IRATE .and. nlangevin .gt. 0) then
        do il=1,nlangevin
           ia=ilangevin(il)

           CALL RANVEL(vx,vy,vz,AMASS(ITY),TEXTL)
cPres   scale the velocity                                                                
           VA(3*ia-2)=colrat*VA(3*ia-2)+sqrt(1-colrat**2)*vx*axinv
           VA(3*ia-1)=colrat*VA(3*ia-1)+sqrt(1-colrat**2)*vy*axinv
           VA(3*ia  )=colrat*VA(3*ia  )+sqrt(1-colrat**2)*vz*axinv

      CALL KINET

      PKNEW=TTOT

      ekdiff=(PKINET-PKNEW)

cd      write(6,*) '#        time = ',time,'   ekdiff = ',ekdiff                           

      RANKIN=RANKIN+(PKINET-PKNEW)
c          this is the change in kinetic energy due to the stoch. coll.                   
      PKINET=PKNEW

      pressure=(pkinet-0.5*virTOT)*ftemp/volm

      if(IPRES) then
cd         volv=0.5*(volv+volapp)+(0.5*stpsz/pmass)*(pressure-pextl)                      
c        Estimate of the time derivative of volume for the next step                      
c        needs to be revised: (just as if the calculation was being started now           
c        using the reassigned velocities as starting velocities):                         
         volapp=volv+stpsz*(pressure-pextl)/pmass

      endif

        enddo

      endif
c                  close `      IF(IRATE .and. nlangevin .gt. 0) '

c             - - - - - - -  - - - - - - - - -  - - - - - - - -                           
c
          if (Iqkmin)  then
cw              write(lunqckmin,*) '   From Mdloop:'
cERB          call QUICKMIN(RAold, VAold, faOld)
              call QUICKMIN
              if (genConstraints) call projectForce
              call rigidatoms(dummy)
          endif
c     zero velocity components perpendicular to the new force
c     for quick minimization, done every timestep.

      call KINET
      pknew = ttot
      rankin = rankin + (pkinet-pknew)
      pkinet = pknew
c       Check if the RMS force has dropped below tolerance, if so jump out of
c       the loop over timesteps:
      if(Iqkmconv) then
         write(lunout,905) tolforce
905      format('#'/'#   ** Convergence reached in Quickmin, ',
     +          ' MAX force < ',g12.4)
         go to 9000
      endif

c                 -------------------------    end quickmin.

      totene=pkinet+ppoten
c             total energy after stochastic collisions.                                   
C   ------------------------------------------------------------------------              
c                                                                                         
c  A new step has been taken, so update the time:                                         
cERB---------------------
cERB Change the way of keeping track of the time to allow a variable time step
cERB      TIME=time0+itime*stpsz
      time = time + stpsz
cERB---------------------
c
cGraphics:
c      call graphics subr. here 

c-----------------------------------------------------------------------
c
C  CALCULATE THE AVERAGES and RMS FLUCTUATIONS: 
                                          
c
c   Add to the interval averaging accumulators: 
                                          
c     First those used to calculate averages over NGRP steps: 
                            
        SUMTEM=SUMTEM+(PKINET-REFKIN)
c     the temperature grid:

cERB Temp. grid ---(start)---
        if (itempgr) then
           write(luntgr, 1328) time,ftemp*sumtem
 1328      format('#'/'#  '/
     +      '#      -----------------------------------------------'/
     +      '#  At time = ',f12.5,' ave temp = ',g15.7)


           azT = 27.d0
           do ix = 1, 3
              do iy = 1, 3
                 do iz = 1, 5
                    avengr(ix, iy, iz) = 0.d0
                    ngr(ix, iy, iz) = 0
                 end do
              end do
           end do
           
           do iAt = 1, nAtms(2)
              ind = 3 * nAtms(1) + 3 * (iAt -1)
              x = ra(ind+1) * ax
              y = ra(ind+2) * ax
              z = ra(ind+3) * ax
              ix = int( (x + ax/2.d0) / (ax/3.d0)) + 1
              iy = int( (y + ay/2.d0) / (ay/3.d0)) + 1
              iz = int( (z + azT/2.d0) / (azT/4.d0)) + 1

              if (iz .le. 5) then
                 ngr(ix, iy, iz) = ngr(ix, iy, iz) + 1
                 ind = 3 * nAtms(1) + 3 * (iAt -1)
                 vel2 = va(ind+1)*va(ind+1) + va(ind+2)*va(ind+2) +
     $                va(ind+3)*va(ind+3) 
                 avengr(ix, iy, iz) = avengr(ix, iy, iz) + 8. * vel2
                 
                 ind = 6 * (iAt -1)
                 vel2 = va(ind+1)*va(ind+1) + va(ind+2)*va(ind+2) +
     $                va(ind+3)*va(ind+3) 
                 avengr(ix, iy, iz) = avengr(ix, iy, iz) + vel2 / 2.d0

                 ind = ind + 3
                 vel2 = va(ind+1)*va(ind+1) + va(ind+2)*va(ind+2) +
     $                va(ind+3)*va(ind+3) 
                 avengr(ix, iy, iz) = avengr(ix, iy, iz) + vel2 / 2.d0
              end if
           end do

           eVtoK = 11604.49987d0 * ax * ax
           do iz = 5, 1, -1
              do iy = 3, 1, -1
                 do ix = 1, 3
                    if (ngr(ix, iy, iz) .ne. 0) then
                       avengr(ix,iy,iz) = avengr(ix,iy,iz) / 
     $                      ngr(ix,iy,iz)
                       avengr(ix,iy,iz) = avengr(ix,iy,iz) * eVtoK
                    end if
                 end do
                 write(luntgr,'(3f15.5)') avengr(1,iy,iz), avengr(2,iy
     $                ,iz), avengr(3,iy,iz)
              end do
              write(luntgr, '()')
           end do

           write(77,*) '--------------------------'
           do iz = 5, 1, -1
              do iy = 3, 1, -1
                 write(77,'(3i10)') ngr(1,iy,iz), ngr(2,iy,iz), ngr(3,iy
     $                ,iz) 
              end do
              write(77, '()')
           end do
        end if
cERB Temp. grid ---(end)---

cERB
cERB        if (itempgr) then
        if (1 .lt. 0) then
cERB
           do 93 ix=1,21
              do 94 iy=1,21
                 do 95 iz=1,21
                    sumtgr(ix,iy,iz)=sumtgr(ix,iy,iz)+tgrid(ix,iy,iz)
                    avengr(ix,iy,iz)=avengr(ix,iy,iz)+natgr(ix,iy,iz)
 95              continue
 94           continue
 93        continue
        end if
        SUMPOT=SUMPOT+(PPOTEN-REFPOT)
        SUMENG=SUMENG+((PKINET+PPOTEN)-REFENG)

c     conserved energy:                                                                   
        ECNSRV=PKINET+PPOTEN+rankin-vdkin+remkin+flwkin

c             add energy exchange to heat bath if stochastic coll. were done.             
c             add energy of atoms that were removed and subtract energy added
c             when new atoms were introduced.
        pistenergy=0.5*pmass*volv**2
     +                                 +pextl*volm
        if(IPRES) ECNSRV=ECNSRV+pistenergy
     +                                    +ranpis
c             add energy of the piston if these are isobaric simulations.                 

        SUMCNS=SUMCNS+(ECNSRV-refeng)
        sumpre=sumpre+(pressure-refpre)
        sumvol=sumvol+(volm-refvol)
c                                                                                         
c      Secondly, add to accumulators used over the whole run:                             
        nsal=nsal+1
        saltem=saltem+pkinet-refkin
        salpot=salpot+ppoten-refpot
        saleng=saleng+totene-refeng
        salcns=salcns+ecnsrv-refeng
        salpre=salpre+pressure-refpre
        salvol=salvol+volm-refvol
c                                                                                         
c      Add to accumulators used to find RMS fluctuations:                                 
        SSQTEM=SSQTEM+(PKINET-REFKIN)**2
        SSQPOT=SSQPOT+(PPOTEN-REFPOT)**2
        SSQENG=SSQENG+(totene-REFENG)**2
        SSQCNS=SSQCNS+(ECNSRV-refeng)**2
        ssqpre=ssqpre+(pressure-refpre)**2
        ssqvol=ssqvol+(volm-refvol)**2
c                                                                                         

c   Calculate averages for the FPI chains:
      if(nFPI .gt. 0) then
c        For each FPI chain (irrespective of what kind of chain it is), get
c          - average position
c          - distribution (standard deviation from average position)
c          - average force
c         all averaged over the images belonging to the same chain.

        do if=1,nFPI
          sumx=0.0
          sumy=0.0   
          sumz=0.0
          sumx2=0.0
          sumy2=0.0
          sumz2=0.0
          sumfx=0.0
          sumfy=0.0
          sumfz=0.0
          indFPI=3*(if-1)
          do im=1,nimFPI(1)
c                assume here that all FPIs have same number of images.
             ind=3*(iFPI(if)+im-2)

cw             write(lunoutFPI,*) '  Loop:  if=',if,'  im=',im
cw             write(lunoutFPI,*) '     Ra = ',RA(ind+1),RA(ind+2),RA(ind+3)

             sumx=sumx+RA(ind+1)*ax
             sumy=sumy+RA(ind+2)*ax
             sumz=sumz+RA(ind+3)*ax
             sumx2=sumx2+(RA(ind+1)*ax)**2
             sumy2=sumy2+(RA(ind+2)*ax)**2
             sumz2=sumz2+(RA(ind+3)*ax)**2
             sumfx=sumfx+FA(ind+1)*ax
             sumfy=sumfy+FA(ind+2)*ax
             sumfz=sumfz+FA(ind+3)*ax
          enddo

c--------------------------------
cERB There are no 300 chains
          if (.FALSE.) then

          if(itag(if) .ge. 300 .and. itag(if) .lt. 399) then
             sumall=abs(sumfx)+abs(sumfy)+abs(sumfz)
             if(sumall .gt. 0.1e-10) then
                write(lunout,*) '#  WARNING from LOOP:  cntr force='
                write(lunout,*) sumfx,sumfy,sumfz
             endif
          endif

          end if
c--------------------------------

          aveposFPI(indFPI+1)=aveposFPI(indFPI+1)+sumx
          aveposFPI(indFPI+2)=aveposFPI(indFPI+2)+sumy
          aveposFPI(indFPI+3)=aveposFPI(indFPI+3)+sumz
          avedistFPI(indFPI+1)=avedistFPI(indFPI+1)+sumx2
          avedistFPI(indFPI+2)=avedistFPI(indFPI+2)+sumy2
          avedistFPI(indFPI+3)=avedistFPI(indFPI+3)+sumz2

          aveforcFPI(indFPI+1)=aveforcFPI(indFPI+1)+sumfx
          aveforcFPI(indFPI+2)=aveforcFPI(indFPI+2)+sumfy
          aveforcFPI(indFPI+3)=aveforcFPI(indFPI+3)+sumfz

cc          aveforcFPI(indFPI+1)=aveforcFPI(indFPI+1)+cntrFA(indFPI+1)
cc          aveforcFPI(indFPI+2)=aveforcFPI(indFPI+2)+cntrFA(indFPI+2)
cc          aveforcFPI(indFPI+3)=aveforcFPI(indFPI+3)+cntrFA(indFPI+3)

        enddo
      endif
c                 close the FPI averaging.


      JGRP=JGRP+1
c           this variable counts how many values have been added to accumultrs.           



c  Check whether group averages should be calculated and printed:                         
      IF(JGRP .eq. NGRP) then
       JGRP=0
c     Now its time to write the averages to the output file:                              
       FACT=1.0/FLOAT(NGRP)
       SUMTEM=(FACT*SUMTEM+REFKIN)/float(natoms-naperm)
c                     correct the temperature for the rigidatoms.                         
c                     atoms by subtracting NAPERM.                                        
c      the temperature grid:                                                              
cERB
cERB        if (itempgr) then
        if (1 .lt. 0) then
cERB
           do 96 ix=1,21
              do 97 iy=1,21
                 do 98 iz=1,21
                    sumtgr(ix,iy,iz)=fact*ftemp*sumtgr(ix,iy,iz)
                    avengr(ix,iy,iz)=fact*avengr(ix,iy,iz)
 98              continue
 97           continue
 96        continue
        end if
c     
      SUMPOT=FACT*SUMPOT+REFPOT
      SUMENG=FACT*SUMENG+REFENG
      SUMCNS=FACT*SUMCNS+refeng
      sumpre=fact*sumpre+refpre
      sumvol=fact*sumvol+refvol
c                                                                                         
C     if(.not. IPRES) then
C        write(lunout,326) time,ftemp*sumtem,sumpre,
C    +                            sumpot,sumeng,sumcns
C     end if

cERB 326  format(f9.3,2g15.7,g15.7,g15.7,g15.7)
C326  format(f9.3,2g15.7,g19.10,g19.10,g15.7)

C Modified by Fer:
C Better looking output
      if( .not. IPRES) then
         write(lunout,FMT='(A,F7.2,2E12.5,3F14.8)')
     &         ' Stp: ',
     &         time,
     &         ftemp*sumtem, sumpre,
     &         sumpot, sumeng, sumcns
      end if

      if(IPRES) write(lunout,360) time,ftemp*sumtem,sumpre,sumvol,
     +                        sumpot,sumeng,sumcns
cERB 360  format(f9.3,g15.7,g15.7,g15.7,g15.7,g15.7,g15.7)
 360  format(f9.3,g15.7,g15.7,g19.7,g19.12,g19.12,g15.7)

cERB      if(ITEMPGR) then
      if(1 .lt. 0) then
        write(luntgr,328) time,ftemp*sumtem
328     format('#'/'#  '/
     +         '#      -----------------------------------------------'/
     +         '#  At time = ',f12.5,' ave temp = ',g15.7)

c    fix tgrid output so as to relax assumptions about the placement of the grid.         
c    program now checks to see if there are no atoms in each slice and then               
c    suppresses output.                                                                   

        do 80 izp = 1, 2*maxgrz+1
          iz = maxgrz - (izp - 1) + 10
c        Check to see if there are any atoms in this layer:                               
          totatgr = 0.0
          do isy = 1, 2 * maxgry + 1
          iy = +maxgry - (isy - 1) + 10
          do isx = 1, 2 * maxgrx + 1
             ix = +maxgrx - (isx - 1) + 10
             totatgr = totatgr + avengr(ix,iy,iz)
          enddo
        enddo
          zgrlow = (iz - 10.5) * tgrlz
          zgrhig = (iz - 9.5) * tgrlz
          write(luntgr,330) iz, zgrlow, zgrhig, totatgr
330       format('#'/
     +           '#   iz = ',i3,'    range:  zmin = ',f7.3,'  zmax = ',
     +           f7.3,'  tot#ofAtoms ',f7.3)
        if(totatgr .lt. 0.1) go to 80

          do 81 iyp=1,2*maxgry+1
           iy=+maxgry-(iyp-1)+10
           ixmin=10-maxgrx
           ixmax=10+maxgrx
           write(luntgr,329) (sumtgr(ix,iy,iz),ix=ixmin,ixmax)
           write(luntgr,344) (avengr(ix,iy,iz),ix=ixmin,ixmax)
329        format(' ',11g9.2)
344        format('#  ',11g9.2)
81         continue
80        continue
c                                                                                         
c        Reset the accumulators:                                                          
         do 290 ix=1,21
          do 291 iy=1,21
            do 292 iz=1,21
              sumtgr(ix,iy,iz)=0.0
              avengr(ix,iy,iz)=0.0
292           continue
291         continue
290        continue
        endif
c                close 'if(ITEMPR)'                                                      


       if(nFPI .gt. 0) then
c    Write out into on FPI into file  outFPI.dat:

        fact=fact/nimFPI(1)
c          assume here that all FPI have same number of images.
cERB        write(lunoutFPI,345) time
345     format('# ',g12.4)
        do if=1,nFPI
          indFPI=3*(if-1)
          aveposFPI(indFPI+1)=aveposFPI(indFPI+1)*fact
          aveposFPI(indFPI+2)=aveposFPI(indFPI+2)*fact
          aveposFPI(indFPI+3)=aveposFPI(indFPI+3)*fact
          sqarg=     fact*avedistFPI(indFPI+1)-aveposFPI(indFPI+1)**2+
     +               fact*avedistFPI(indFPI+2)-aveposFPI(indFPI+2)**2+
     +               fact*avedistFPI(indFPI+3)-aveposFPI(indFPI+3)**2
          if(sqarg .gt. 0) then
            stdev=sqrt(sqarg)
          else
            stdev=-99999.9
cERB            write(lunoutFPI,*)' Warning: sqarg in calc stdev < 0,',sqarg
          endif
          aveforcFPI(indFPI+1)=aveforcFPI(indFPI+1)*fact
          aveforcFPI(indFPI+2)=aveforcFPI(indFPI+2)*fact
          aveforcFPI(indFPI+3)=aveforcFPI(indFPI+3)*fact

cERB          write(lunoutFPI,346) if,(aveposFPI(indFPI+k),k=1,3),stdev,
cERB     +                  (aveforcFPI(indFPI+k),k=1,3),sumeng
346       format('   FPI',i3,'   ',3g12.4,g16.4,3g12.4,g16.4)
        enddo

c    Zero accumulators:
        do if=1,MAXFPI
         ind=3*(if-1)
         aveposFPI(ind+1)=0.0
         aveposFPI(ind+2)=0.0
         aveposFPI(ind+3)=0.0
         avedistFPI(ind+1)=0.0
         avedistFPI(ind+2)=0.0
         avedistFPI(ind+3)=0.0
         aveforcFPI(ind+1)=0.0
         aveforcFPI(ind+2)=0.0
         aveforcFPI(ind+3)=0.0
        enddo
       endif
c           close   'if(nFPI .gt. 0) then '

        SUMTEM=0.0
        SUMPOT=0.0
        SUMENG=0.0
        SUMCNS=0.0
        sumpre=0.0
        sumvol=0.0

      endif
c           close 'if(JGRP.eq.NGRP)'.                                                     

c                           ----------------------------                                                                               
c   The average coordinates:                                                              
         javecoo=javecoo+1
       if(javecoo .eq. nstaveRA) then
            javecoo=0
          fact=1./float(navecoo)
c                                                                                         
cq        write(lunout,501) navecoo                                                       
cq501     format('#    navecoo = ',i10,                                                    
cq     +         '   Ra(1), sumRA(1), aveRA(1) = ',3g15.6)                                
c	                                                                                        
          do 293 i=1,3*natoms
             aveRA(i)=fact*sumRA(i)
             sumRA(i)=0.0
293            continue
c         PBC: Make sure the average coordinate lands within the simulation box:          
c            Rectangular simulation cell:                                                 
             DO I=1,NATOMS
                ind=3*(i-1)
c               first the y coordinate:                                                   
                  IF(aveRA(ind+2) .GT. ayhalfsc) then
                     aveRA(ind+2)=aveRA(ind+2)-aysc
                  endif
                  IF(aveRA(ind+2) .LE. -ayhalfsc) then
                     aveRA(ind+2)=aveRA(ind+2)+aysc
                  endif
c               then the x coordinate:                                                    
                  IF(aveRA(ind+1) .GT. axhalfsc) then
                     aveRA(ind+1)=aveRA(ind+1)-axsc
                  endif
                  IF(aveRA(ind+1) .LE. -axhalfsc) then
                     aveRA(ind+1)=aveRA(ind+1)+axsc
                  endif
c               finally the z coordinate:                                                 
                  IF(aveRA(ind+3) .GT. azhalfsc) then
                     aveRA(ind+3)=aveRA(ind+3)-azsc
                  endif
                  IF(aveRA(ind+3) .LE. -azhalfsc) then
                     aveRA(ind+3)=aveRA(ind+3)+azsc
                  endif
             enddo
c                                                                                         
          navecoo=0
c            the array aveRA now contains the most uptodate average coordnts.	            
          endif
c                                                                                         
c   -----------------------------------------------------------------                     
c                                                                                         
c  Dump an Intermediate configuration:                                                    
      if(IDUMP) then
	 
        JDUMP=JDUMP+1
        if(JDUMP .ge. NWAIT) IWAIT=.false.
c                                   then no more waiting.                                 
        if( .not. IWAIT .and. JDUMP .ge. NDUMP) then
c                                                                                         
           call conout(luncm,indvel,ax)
           JDUMP=0
        endif
      endif
c                                                                                         
c   -----------------------------------------------------------------                     
c                                                                                         
c  Write info to the viscos.dat file:                                                     
      if(IVISC) then
        JVISC=JVISC+1
        if(JVISC .ge. NWVISC) IWVISC=.false.
c                                   then no more waiting.                                 
        if( .not. IWVISC .and. JVISC .ge. NVISC) then
c                                                                                         
           NVISCPT=NVISCPT+1
c        find energy per atom:                                                            
         enperat=totene/float(NATOMS)
           temper=ftemp*pkinet/float(natoms-naperm)

           write(lunvis,480) NVISCPT,time,temper,pressure,enperat
480        format('#',i5,f10.4,g22.12,g22.12,g22.12)

           sumxpx=0.
           sumypy=0.
           sumzpz=0.
           sumxde=0.
           sumyde=0.
           sumzde=0.
           sumxpy=0.
           sumypx=0.
           sumxpz=0.
           sumzpx=0.
           sumypz=0.
           sumzpy=0.

c        Calculate the sums over all atoms:                                               
c        First loop over components:                                                      
           cooscale=ax
           isumat=0
           do 48 itype=1,NATYPE
             am=amass(itype)
             do 49 k=isumat+1,isumat+natms(itype)
c                      use coordinates that have not been subject to PBC:                 
                 x=raNoPBC(3*K-2)*cooscale
                 y=raNoPBC(3*K-1)*cooscale
                 z=raNoPBC(3*K)*cooscale
c                     find momentum for each atom:                                        
                 Px=am*VA(3*K-2)*cooscale
                 Py=am*VA(3*K-1)*cooscale
                 Pz=am*VA(3*K)*cooscale
c              find deviation in total energy from average total energy for each:   
                 De=(0.5*(Px**2+Py**2+Pz**2)/am+Potpat(k))-enperat

             sumxpx=sumxpx+x*Px
             sumypy=sumypy+y*Py
             sumzpz=sumzpz+z*Pz

             sumxde=sumxde+x*De
             sumyde=sumyde+y*De
             sumzde=sumzde+z*De

             sumxpy=sumxpy+x*Py
             sumypx=sumypx+y*Px
             sumxpz=sumxpz+x*Pz
             sumzpx=sumzpx+z*Px
             sumypz=sumypz+y*Pz
             sumzpy=sumzpy+z*Py

cq                     write(lunvis,485) k,x,y,z                                          
cq485                  format('# Atom  ',i4,':   x,y,z = ',3g22.12)                        
cq                     write(lunvis,483) x*px,y*py,z*pz,x*de,y*de,z*de                    
cq483                  format('#      ',3g22.14/'      ',3g22.14)	        	                
cq                     write(lunvis,484) x*py,y*px,x*pz,z*px,y*pz,z*py                    
cq484                  format('#      ',3g22.14/'      ',3g22.14)	        	                

49               continue
            isumat=isumat+natms(itype)
48            continue

c          divide by the volume:                                                          
            sumxpx=sumxpx/volm
            sumypy=sumypy/volm
            sumzpz=sumzpz/volm
            sumxde=sumxde/volm
            sumyde=sumyde/volm
            sumzde=sumzde/volm
            sumxpy=sumxpy/volm
            sumypx=sumypx/volm
            sumxpz=sumxpz/volm
            sumzpx=sumzpx/volm
            sumypz=sumypz/volm
            sumzpy=sumzpy/volm

           write(lunvis,481) sumxpx,sumypy,sumzpz,sumxde,sumyde,sumzde
481        format('        ',3g24.14/'        ',3g24.14)
           write(lunvis,482) sumxpy,sumypx,sumxpz,sumzpx,sumypz,sumzpy
482        format('        ',3g24.14/'        ',3g24.14)

           JVISC=0
        endif
      endif
c                                                                                         
c   -----------------------------------------------------------------                     
C                                                                                         
C   WRITE INSTANTANEOUS VALUES:                                                           
C                                                                                         
      IF(INSTA) then
        JINST=JINST+1
        IF(JINST .ge. NINST) then
c                                                                                         
         temper=ftemp*pkinet/float(natoms-naperm)
c                                                                                         
         if(IPRES) write(lunins,321) time,temper,pressure,volm,
     +                                     ppoten,totene,ecnsrv
321      format(' ',7g22.14)
         if(.not.IPRES) write(lunins,327) time,
cc     +                   ra(1)*ax,ra(2)*ax,ra(3)*ax,
     +                   temper,pressure,ppoten,totene,ecnsrv                   
327      format(' ',10g22.14)
c                                                                                         
         JINST=0
       endif
      endif
c   -----------------------------------------------------------------                     

C   WRITE trajectory of a selected atom:                                                           

      IF(Iattraj) then
        Jstatwr=Jstatwr+1
        IF(Jstatwr .ge. Nstatwr) then
         ind=3*(indatwr-1)
c Comment all this part, and uncomment commented lines to track top layer
         write(lunattr,232) time,
     +         ra(ind+1)*ax,ra(ind+2)*ax,ra(ind+3)*ax,
     +         va(ind+1)*ax,va(ind+2)*ax,va(ind+3)*ax,
     +         fa(ind+1)*ax,fa(ind+2)*ax,fa(ind+3)*ax,
cc     +         potenpat(indatwr),
     +         ppoten
 232     format(' ',f10.5,'  ',3f9.4,'  ',3f9.4,'  ',3f10.5,2g22.12)

c Line added to trace the potential of the surface
         write(99, *) ra(ind+3)*ax, -ppoten
c ---------------------------------------------ERB

c   
cERB Modification to keep track of the top layer
cERB         write(lunattr, '(f10.5,13f8.5)') time,
cERB     $        (ra(3*(170+i)), i = 1, 12), ra(3*255)
cERB
         Jstatwr=0

C Added by Fer.
C Write the number of molecules
C       write(87,FMT='(I5)') 3*nMolecules

C Write the title
C       write(87,FMT='(A)') ' '

C Write the XYZ coordinates
C       do i=1,nMolecules

C Get the coordinates for molecule i
C         Pt_i1 = 6*(i-1)
C         Pt_i  = 3*(i-1+2*nMolecules)
C         Pt_i2 = 3+6*(i-1)

C         do p=1,3
C           Ri1(p) = box(p)*ra(Pt_i1 + p)
C           Ri(p)  = box(p)*ra(Pt_i  + p)
C           Ri2(p) = box(p)*ra(Pt_i2 + p)
C           Ri1(p) = ax*ra(Pt_i1 + p)
C           Ri(p)  = ax*ra(Pt_i  + p)
C           Ri2(p) = ax*ra(Pt_i2 + p)
C         end do

C         write(87,FMT='(A,3F16.10)') '  H',(Ri1(p),p=1,3)
C         write(87,FMT='(A,3F16.10)') '  O',(Ri(p),p=1,3)
C         write(87,FMT='(A,3F16.10)') '  H',(Ri2(p),p=1,3)

C       end do

C Added by Fer.
C Dump important stuff
        write(89) time, uTot, (ax*Ra(i),ax*Va(i),i=1,9*nMolecules)

       endif
      endif

c   --------------------------------------------------------------------------
c      ---------------------------------------------------------------------

cSurface options:
c
c  REMOVE ATOMS FLYING AWAY:
c     Check if an atom has moved far from the surface, farther than -ZFLAWY
c     and is heading away from it. If such an atom is found, if is removed.
         zflawy=zsurf+1.1*rskinmax
c                   an atom will be removed if it is more than 10% above the
c                   highest atom (zsurf) and is moving away.
         nrmnow=0
         narnow=0
         do 118 i=1,natoms
            ind=3*(i-1)
            if(ra(ind+3) .gt. zflawy 
     +                          .and. va(ind+3) .gt. 0.1e-06) then
              nflawy=nflawy+1
              nrmnow=nrmnow+1
c                   number of atoms that will be removed now.
              ityp=1
              if(i .gt. natms(1)) ityp=2
	      if(ityp .eq. 1) narnow=narnow+1
c                   counts the number of A atoms that are removed.
              iatfla(nflawy+nrmnow)=i
c                   store the index of the atoms that are removed.
c
c            Kinetic energy of this atom is lost, correct ECONST through RANKIN:
              delv2=0.5*amass(ityp)*(va(ind+1)**2+
     +                 va(ind+2)**2+va(ind+3)**2)
              flwkin=flwkin+delv2
c
              write(lunout,500) i,time,(ra(ind+k),k=1,3),zflawy,
     +                      (va(ind+kk),kk=1,3),delv2,zbott,zmaxperm
500           format(/'   *** Atom ',i4,' is flying away,',
     +              '    time=',f10.3/
     +        '    x,y =',2g12.4,'  z =',g12.4,' > zflawy =',g12.4/
     +              '    Vx, Vy, Vz =',3g12.4,
     +              '  delV2 =',g12.4/
     +              '   zbott = ',g12.4,'  zmaxperm = ',g12.4)
c
cG77              call removeatom(i)
c
            endif
118         continue
c
c
c             ---------------------------------------------------
c
c   ATOM DEPOSITION:
c     Check whether an atom should be deposited on the surface:
        jvdep=jvdep+1
c           this variable counts the number of timesteps since last ATOM DEP.
      if(Ideposit .and. jvdep .eq. ntotvd) then
c            send a new atom towards the surface:
          jvdep=0
c            reset counter.
c
cd              write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))
c
c      Just before adding a new atom the position of the surface boundary is
c      updated (so that T grid is reasonably well positioned). Find MAX(z):
          zlarge=zbott
          do 191 i=1,natoms
            if(ra(3*i) .gt. zlarge) then
	      zlarge=ra(3*i)
	      indzsurf=i
	    endif
191         continue
          zsurf=zlarge
c
c    Decide which kind of atom to deposit:
        itypvd=1
c            1 means that the atom will be an A atom.
        if(ranb(iseed1,iseed2) .gt. afract) itypvd=2
c
c    Decide what the (x,y) position of the new atom will be: (assume rectangle)
       scalern=axsc
       if(aysc .gt. axsc) scalern=aysc
c           make sure the square of radom impacts covers the simulation box.
c       Now find a pair of random numbers which give impact parameter inside box
c       Maximum of 1000 tries allowed
        do 660 irantry=1,1000
          x0vd=(ranb(iseed1,iseed2)-0.5)*scalern
          y0vd=(ranb(iseed1,iseed2)-0.5)*scalern
          if(x0vd .gt. axhalfsc) go to 660
          if(x0vd .le. -axhalfsc) go to 660
          if(y0vd .gt. ayhalfsc) go to 660
          if(y0vd .le. -ayhalfsc) go to 660
c           here an impact parameter has been found that lies inside.
          go to 661
660       continue

        write(lunout,*) '# WARNING Loop:  cant find impact param inside'
661     continue

cd        write(lunout,595) zsurf,x0vd,y0vd
cd595     format('    zsurf= ',g20.10,'   x0vd,y0vd= ',2g15.6)

c    Decide what the z position of the new atom will be:
c        Find the highest atom within a cylinder of radius RSKIN around the
c        (x,y) position and start the new atom a  DZ0VD above that 
c        atom:
c          (Use RSKIN rather than RCUT just to be safe (vibrations)).
c          Simply use the larger of the two possible values of RSKIN:
            if(itypvd .eq. 1) rskinvdmax2=dmax1(rskin2(1),rskin2(3))
            if(itypvd .eq. 2) rskinvdmax2=dmax1(rskin2(2),rskin2(3))
c
cd          write(lunout,597) rskinvdmax2
cd597       format('    rskinvdmax2 = ',g20.10)

	    zmaxcyl=zmaxperm
            do 401 i=1,natoms
	      if(zmaxcyl .gt. ra(3*i)) go to 401
	        delx=ra(3*(i-1)+1)-x0vd
	        dely=ra(3*(i-1)+2)-y0vd
c
c          PBC:  Periodic boundary conditions:
c              First the y coordinate:
                 IF(dely .GT. ayhalfsc) dely=dely-aysc
                 IF(dely .LT.-ayhalfsc) dely=dely+aysc
c             Then the x coordinate:
                 IF(delx .GT. axhalfsc)delx=delx-axsc
                 IF(delx .LT.-axhalfsc) delx=delx+axsc
c
		xydist2=delx*delx+dely*dely
		if(xydist2 .gt. rskinvdmax2) go to 401
		  zmaxcyl=ra(3*i)
401             continue		
c
            z0vd=dz0vd+zmaxcyl
c
            if((natoms+1) .gt. maxatoms) then
                 write(lunout,*) '  Too many atoms, natoms = ',natoms,
     +                           '  maxatoms = ',maxatoms
                 go to 100

            endif
c
            call addatom(itypvd,x0vd,y0vd,z0vd,0.,0.,v0vd,navecoo)
c
            natdep=natdep+1
c
            vdkin=vdkin+0.5*amass(itypvd)*v0vd**2
c                      keep track of the energy that is added to the
c                      system by the VD atoms and add that to ECNSRV
c
             write(lunout,590) itypvd,natoms,natms(1),time,zmaxcyl*ax,
     +                     x0vd*ax,y0vd*ax,(dz0vd+zmaxcyl)*ax,v0vd*ax
590          format('#'/
     +              '# Deposit type ',i1,' atom:  Ntot =',i5,' Na ='
     +              ,i5,' time =',f9.3,' zmaxcyl',g12.4/
     +              '    x0 =',g12.4,' y0 =',g12.4,' z0 =',g12.4,
     +              '   V0z =',g12.4)
c
cd              write(lunout,*)'  After atom deposition:  naperm=',naperm	  
cd              write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))
c                 ------------------------------
cc        go to 900
c
c    Lift the PERMAFROST up by one atom:
c        Find the lowest that is not in the permafrost:
            indaddper=0
	    znextperm=zsurf+0.0001
            do 196 i=1,natoms
               do k=1,naperm
                 if(i .eq. iperm(k)) go to 196
c                          limit this to atoms that are NOT in permafrost.	
               enddo
c            It is important here to use the AVERAGE coordinate so that
c            atoms will not get trapped in permafrost because of vibration:
               if(avera(3*i) .ge. znextperm) go to 196
		 znextperm=avera(3*i)
		 indaddper=i
c
cd                write(lunout,594) i,k,naperm,znextperm
cd594             format('         i,k,naperm,znextperm = ',3i6,g15.6)
c	   	   
196            continue

        go to 900
c              -----------------------
c         add this atom to the list of permafrost atoms:
            naperm=naperm+1
            naddper=naddper+1
	    iperm(naperm)=indaddper
            itag(indaddper)=1
	    itypaddper=1
	    if(indaddper .gt. natms(1))itypaddper=2
c         update zmaxperm:
            zmaxperm=znextperm
c         keep track of the displacement as the atoms are moved to average coo:
            disp2=(ra(3*(indaddper-1)+1)-aveRA(3*(indaddper-1)+1))**2+
     +            (ra(3*(indaddper-1)+2)-aveRA(3*(indaddper-1)+2))**2+
     +            (ra(3*(indaddper-1)+3)-aveRA(3*(indaddper-1)+3))**2
            disp=sqrt(disp2)
            if(disp .gt. dispmaxaddper) dispmaxaddper=disp
            disp2addper=disp2addper+disp2
c
            ra(3*(indaddper-1)+1)=aveRA(3*(indaddper-1)+1)
            ra(3*(indaddper-1)+2)=aveRA(3*(indaddper-1)+2)
            ra(3*(indaddper-1)+3)=aveRA(3*(indaddper-1)+3)
c                  place this atom at its AVERAGE coordinates in the permafrost
c
c        Kinetic energy of this atom is lost, later correct ECONST. (see 
c                 beginning of dynamics loop)
c
c        Now update 'zmindynat' which is the z coordinate of the lowest 
c        dynamical atom and update 'ndynbelper' which is the number of dynamical
c        atoms below the highest permafrost atom:       
            ndynbelper=0
            indlowdyn=0
	    zcoo=zsurf+0.0001
            do 197 i=1,natoms
               do k=1,naperm
                 if(i .eq. iperm(k)) go to 197
c                          limit this to atoms that are NOT in permafrost.	
               enddo
c           check the AVERAGE coordinate of this dynamical atom:
               if(avera(3*i) .lt. zmaxperm) ndynbelper=ndynbelper+1 
               if(avera(3*i) .ge. zcoo) go to 197
		 zcoo=avera(3*i)
		 indlowdyn=i
c
cd                write(lunout,594) i,k,naperm,zcoo
c	   	   
197            continue
             zmindynat=zcoo
c
          write(lunout,591) itypaddper,indaddper,zmaxperm,naperm,
     +             zmindynat,ndynbelper
591       format(' Add atom of type ',i2,' to permafrost: indaddper = ',i8,
     +           ' zmaxperm = ',g15.5/'    naperm = ',i5,
     +           '  zmindynat =',g12.4,' ndynbelper',i5)
c
cd         write(lunout,599) VA(3*indaddper-2),VA(3*indaddper-1),VA(3*indaddper) 
cd599      format('     Velocity components:   x,y,z= ',3g15.6)

c                  ---------------------------------
c
c          go to 900		    
c    Remove the lowest atom if the PERMAFROST is thick enough (at least as
c         thick as the SKINdiameter):
c         The lowest atom is atom 'indzbott' at ZBOTT
c         The lowest dynamical atom is at  'zmindynat'
            if((zmindynat-zbott) .gt. rskinmax) then
	        ninfraperm=ninfraperm+1
c               rinfraperm(3*(ninfraperm-1)+1)=ra(3*(indzbott-1)+1)
		rinfraperm(3+1)=ra(3*(1)+1)
		rinfraperm(3*(ninfraperm-1)+2)=ra(3*(indzbott-1)+2)
		rinfraperm(3*(ninfraperm-1)+3)=ra(3*(indzbott-1)+3)
c
                itype=1
                if(indzbott .gt. natms(1)) itype=2
                itinfraperm(ninfraperm)=itype
c
cG77                call removeatom(indzbott)
c	  
                 write(lunout,593) indzbott,zbott,ninfraperm,naperm
593              format(' Remove atom ',i4,' from perma:',
     +                ' z = ',f7.4,' ninfraperm = ',i4,' naperm= ',i4)
c
c       Find where the lowest atom in permafrost is now: (renew ZBOTT)
             zsmall=zsurf
             do i=1,naperm
	         if(ra(3*iperm(i)) .lt. zsmall) then
		    zsmall=ra(3*iperm(i))
		    indzbott=iperm(i)
		 endif
	      enddo
	      zbott=zsmall
c	      
           write(lunout,596) natoms,natms(1),natms(2),zbott,indzbott,naperm
596        format('    natoms,#A,#B=',3i5,' nzbott=',
     +                  g12.4,' atom =',i5,' naperm=',i5/)
c
            endif
c                  close  'if((znearperm-zbott) .gt. rskinmax)' 
900      continue
c	     
         endif
c                   close 'if(Ideposit .and. jvdep .eq. ntotvd)'.
c
c end Surface options.

c   ------------------------------------------------------------------                                       
c    Before taking next step:                                                             
c     increment the bath parameters if desired: (linear change in time of T or P)        
        if(ILINCH) then
          if(indpara .eq. 1) TEXTL=TEXTL+ratlinch*stpsz
          if(indpara .eq. 2) PEXTL=PEXTL+ratlinch*stpsz
        endif

100   continue

      call saddlefind

c$$$      open (76, file = 'rms.dat', status = 'unknown')
c$$$      do iatom = nHydrogens+1, nHydrogens+nOxygens
c$$$         zAve = zCoord(iatom-nHydrogens)/ntotstps * ax
c$$$         z2Ave = zCoord2(iatom-nHydrogens)/ntotstps * ax * ax
c$$$         rms = sqrt(z2Ave - zAve**2) 
c$$$         write (76, '(i5, 3f15.8)') iatom, rms, zAve, z2Ave
c$$$      end do
c$$$      close (76)
c$$$
c$$$c     Save the histogram
c$$$      open (76, file = 'histogram.dat', status = 'unknown')
c$$$      write(76,'(A)') '# Density histogram '
c$$$      write(76,'(A)') '# bin number    z       rho '
c$$$      do i = 1, 1000
c$$$         z = z0 + (i-0.5) * dz
c$$$         write(76, '(i6, f12.5, i15)') i, z, histogram(i)
c$$$      end do
c$$$      close(76)

c    -------------------------------------------------------------------                  
c                        Close Loop over timesteps.                                       
c    -------------------------------------------------------------------                  
c    -------------------------------------------------------------------                  
c    -------------------------------------------------------------------
c                                                                                         
9000  CONTINUE                                                                          
C                                                                                         
c                                                                                         
c  Write properties of final configuration to output file:                                
      tempfi=ftemp*pkinet/float(natoms-naperm)
      if(.not.IPRES) 
     +   write(lunout,235) time,tempfi,pressure,ppoten,totene,ecnsrv
      if(IPRES) 
     +  write(lunout,235) time,tempfi,pressure,volm,ppoten,totene,ecnsrv
cERB 235  format('#'/'# Finally:  '/'#',f8.3,6g15.7)
 235  format('#'/'# Finally:  '/'#',f9.3,2g15.7, 2g19.10, 2g15.7)

c                                                                                         
      tv2perm=0.5*v2perm*ftemp/float(natoms-naperm)
c                  convert the sum over V**2 into temperature.                            
cd      write(lunout,234) pkinet,v2perm,tv2perm                                           
cd234   format('#   pkinet = ',g20.10,' v2perm = ',g15.6,' tv2 =',                         
cd     +         g15.6)                                                                   
c                                                                                         
cd      write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))                             
cd233   format('#   z component of velocities for type',i1,'atoms:'/                       
cd     +        '#',5g15.6/'#',5g15.6/'#',5g15.6/'#',5g15.6)                                              

c   Find RMS fluctuations for the whole run and print out:                                
          fact=1./float(nsal)
          saltem=saltem*fact
          salpot=salpot*fact
          saleng=saleng*fact
          salcns=salcns*fact
          salpre=salpre*fact
          salvol=salvol*fact
c                                                                                         
          SSQTEM=FACT*SSQTEM-saltem**2
          SSQTEM=ftemp*SQRT(ABS(SSQTEM))/float(NATOMS-naperm)
          SSQPOT=FACT*SSQPOT-(salPOT)**2
          SSQPOT=SQRT(ABS(SSQPOT))
          SSQENG=FACT*SSQENG-(salENG)**2
          SSQENG=SQRT(ABS(SSQENG))
          SSQCNS=FACT*SSQCNS-(salCNS)**2
          SSQCNS=SQRT(ABS(SSQCNS))
          ssqpre=fact*ssqpre-salpre**2
          ssqpre=sqrt(abs(ssqpre))
          ssqvol=fact*ssqvol-salvol**2
          ssqvol=sqrt(abs(ssqvol))
c                                                                                         
      avetem=ftemp*(saltem+refkin)/float(natoms-naperm)
      avepot=salpot+refpot
      aveeng=saleng+refeng
      avecns=salcns+refeng
      avepre=salpre+refpre
      avevol=salvol+refvol
c                                                                                         
      if(IPRES) 
     +    write(lunout,244) avetem,avepre,avevol,avepot,aveeng,avecns,
     +             ssqtem,ssqpre,ssqvol,ssqpot,ssqeng,ssqcns
244   format('#'/'# Whole run:',
     +        '    Temprtr:      Pressure:       Volume:',
     +        '     Potntl.En:',
     +        '    Tot.Energy:       Econst:   '/
     +        '#  Average:',1x,g12.3,4g15.7,g14.6/
     +        '#  RMS fluct:',g12.3,4g15.3,g12.2)
      if(.not. IPRES) 
     +    write(lunout,241) avetem,avepre,avepot,aveeng,avecns,
     +             ssqtem,ssqpre,ssqpot,ssqeng,ssqcns
241   format('#'/'# Whole run:',
     +        '   Temprtr:      Pressure:     Potntl.En:',
     +        '    Tot.Energy:       Econst:'/
     +        '#  Average:',1x,4g15.7,g14.6/
     +        '#  RMS fluct:',g12.3,3g15.3,g12.2)
c                                                                                         
c                                                                                         
c       calculate the average number of steps between updates of neighbor list:           
      if(icntud .gt. 0 .and. icntud .lt. MAXNUPTDL) then
         aveupd=0.0
         do iu=1,icntud
           aveupd=aveupd+nuptdl(iu)
         enddo
           aveupd=aveupd/float(icntud)
      else
         aveupd=nstupd
      endif
      write(lunout,245) aveupd
245     format('#'/'#  Average number of steps between updates of',
     +          ' neighborlist=',g12.4)

cq      write(lunout,243)                                                                 
cq243   format('#'/'#  Number of steps between updates:')                                     
cq      write(lunout,242) (nuptdl(i),i=1,icntud)                                          
cq242   format('# ',                                                                       
cq     +        20i4/'#',20i4/'#',20i4/'#',20i4/
c      +        '#',20i4/'#',20i4/'#',20i4/'#',20i4/'#',20i4)                             
c                                                                                         
cq      write(lunout,262) maxneb,maxnpr                                                   
cq262   format('#'/'#  At the end of loop:   maxneb, maxnpr = ',2i15)                         
c                                                                                         
c                                                                                         
c     -------------------------------------------------------------                       
c    Write the final configuration:                                                       
cAllinCh:
         write(6,*) ' Write the final configuration:'
        if(indallinchains .eq. 1) then
           do ia=1,nimFPI(1)
             itag(ia)=ifirsttag
           enddo
        endif

       open(lunco,file='co.con')
cERB       if(IQKMIN) then
cERB           CALL CONOUT(lunco,-1,ax)
cERB       else
           CALL CONOUT(lunco,1,ax)
cERB       endif
c                                                                                         
       if(IAVECOO) then
c        Write the final AVErage configuration and add the atoms removed from             
c        permafrost (if the total number of atoms is greater than maxatoms the            
c        VA and eventually the FA arrays will be overwritten, but thats OK now.           
c        First component one:                                                             
         do i=1,3*natms(1)
           RA(i)=aveRA(i)
           enddo
         ninfraA=0
         do i=1,nAinfraperm
             RA(3*(ninfraA+natms(1)-1)+1)=rinfraperm(3*(i-1)+1)
             RA(3*(ninfraA+natms(1)-1)+2)=rinfraperm(3*(i-1)+2)
             RA(3*(ninfraA+natms(1)-1)+3)=rinfraperm(3*(i-1)+3)
         enddo
c        Total number of component 1 atoms:                                               
         na=natms(1)+nAinfraperm
c        Then component two:                                                              
         do i=1,natms(2)
           RA(3*(na+i-1)+1)=aveRA(3*(i+natms(1)-1)+1)
           RA(3*(na+i-1)+2)=aveRA(3*(i+natms(1)-1)+2)
           RA(3*(na+i-1)+3)=aveRA(3*(i+natms(1)-1)+3)
         enddo
         ninfraB=ninfraperm-nAinfraperm
         ninfraA=nAinfraperm
         ind=natoms+nAinfraperm
         do i=1,ninfraB
             RA(3*(i+ind-1)+1)=rinfraperm(3*(i+ninfraA-1)+1)
             RA(3*(i+ind-1)+2)=rinfraperm(3*(i+ninfraA-1)+2)
             RA(3*(i+ind-1)+3)=rinfraperm(3*(i+ninfraA-1)+3)
             itag(i+ind)=3
          enddo
c         Total number of component 1 atoms:                                              
          nb=natms(2)+ninfraB
          if(ninfraperm .ne. (ninfraA+ninfraB)) write(lunout,*)
     +            '#  WARNING:   ninfraperm .ne. (ninfraA+ninfraB)'
c                                                                                         
c                                                                                         
          natms(1)=na
          natms(2)=nb
          natoms=na+nb
          naperm=naperm+ninfraA+ninfraB
          CALL CONOUT(lunaco,-1,ax)
c                                                                                         
        endif
c                     close 'if(IAVECOO) then'
c     -------------------------------------------------------------        
c  --------------------------------------------------------------------------

cFPI:
c  Get energy for each chain image separately and find the displacement between
c  images in each chain,
c  write out the coordinates of the highest energy image into file `cmaxim.con'
c  if nudging has been turned on:

      if(nFPI .gt. 0) then
        enmax=-999999.99

        write(lunout,601)
cERB 601     format('#'/'#    Energy of FPI images:'/'#    image #: ',
cERB     +      '        pot. energy:  relative pot:  kinetic en:  ',
cERB     +      'RMSforce: MAXforce:')
601     format('#'/'#    Energy of FPI images:'/'#    image #: ',
     +      '  pot. energy:     relative pot:  kinetic en:  ',
     +      'RMSforce: Path Coor')

cERB        write(lunoutFPI,602) nFPI
602     format('#'/'#    Displacement of FPI images along the chain:'/
     +      '#    image #:    FPI chain 1 ... ',i4,'   Potenergy: ')

        nFPIold=nFPI
        natinch=nFPI*nimFPI(1)
        natnotinch=natoms-natinch

        natomsold=natoms
        natoms=natomsold-nFPI*(nimFPI(1)-1)
        natinch=nFPI*nimFPI(1)
        do it=1,MAXCOMP
          natmsold(it)=natms(it)
          natms(it)=natms(it)-nFPIcom(it)*(nimFPI(1)-1)
        enddo

c     copy all coordinates and velocities into working arrays, rinfraperm and           
c     rdispl:                                                                            
        do ia=1,natomsold
          ind=3*(ia-1)
          rinfraperm(ind+1)=ra(ind+1)
          rinfraperm(ind+2)=ra(ind+2)
          rinfraperm(ind+3)=ra(ind+3)
          rdispl(ind+1)=va(ind+1)
          rdispl(ind+2)=va(ind+2)
          rdispl(ind+3)=va(ind+3)
          itagst(ia)=itag(ia)
        enddo
 
cERB          write(lunoutFPI,*)'  image#:  displacement:       ... difpot'
cERB          write(lunoutFPI,*)'           Chain1:   Chain2:  Chain3: ...'  

c            --------

cAllinChain Transitory "solution". The output doesn't work, so skip the
c        energy of each image.        
cAllinChain
        if (indallinchains .eq. 0) then 
cAllinChain
        do 600 image=1,nimFPI(1)
c              asuming here that all FPI have same number of images.                      
          itot=0
          noldFPIat=0
          do ia=1,natomsold
            if(itagst(ia) .lt. 100) then
              itot=itot+1
              ind=3*(itot-1)
              ra(ind+1)=rinfraperm(3*(ia-1)+1)
              ra(ind+2)=rinfraperm(3*(ia-1)+2)
              ra(ind+3)=rinfraperm(3*(ia-1)+3)
              itag(itot)=itagst(ia)

              if(Iqkmin) then
                 va(ind+1)=0.0
                 va(ind+2)=0.0
                 va(ind+3)=0.0
              else                 
                 va(ind+1)=rdispl(3*(ia-1)+1)
                 va(ind+2)=rdispl(3*(ia-1)+2)
                 va(ind+3)=rdispl(3*(ia-1)+3)
              endif

            else

c            this atom corresponds to a FPI image in some chain, check whether            
c            it is the right image:                                                       
              ip=numFPI(itagst(ia))
              imia=ia-iFPI(ip)+1
              if(image .eq. imia) then
c          include this atom in the config:
                itot=itot+1
                ind=3*(itot-1)
                ra(ind+1)=rinfraperm(3*(ia-1)+1)
                ra(ind+2)=rinfraperm(3*(ia-1)+2)
                ra(ind+3)=rinfraperm(3*(ia-1)+3)

               itag(itot)=0
               if(image .eq. 1 .or. image .eq. nimFPI(1))itag(itot)=1
               ioldFPIat(noldFPIat)=itot
               noldFPIat=noldFPIat+1

c          Find the displacement of this image from the previous one:
                if(imia .eq. 1) then 
                  do iff=1,nFPIold
                    displim(ip)=0.0
                  enddo
                else
                  displim(ip)=displim(ip)+
     +                      sqrt((FPIprevim(3*(ip-1)+1)-ra(ind+1))**2+
     +                           (FPIprevim(3*(ip-1)+2)-ra(ind+2))**2+
     +                           (FPIprevim(3*(ip-1)+3)-ra(ind+3))**2)
                endif
c                     close `if(imia .eq. 1)  ...'

                FPIprevim(3*(ip-1)+1)=ra(ind+1)
                FPIprevim(3*(ip-1)+2)=ra(ind+2)
                FPIprevim(3*(ip-1)+3)=ra(ind+3)

                if(Iqkmin) then
                 va(ind+1)=0.0
                 va(ind+2)=0.0
                 va(ind+3)=0.0
                else                 
                 va(ind+1)=rdispl(3*(ia-1)+1)
                 va(ind+2)=rdispl(3*(ia-1)+2)
                 va(ind+3)=rdispl(3*(ia-1)+3)
                endif

              endif
c                  close `if(image .eq. imia) ... '

            endif
c                close `if(itagst(ia) .lt. 100) ...'

          enddo
c              close `do ia=1,natomsold ... '
c           -------------

c       Erase all indication of FPI atoms:  (comment out: done already now)
cc          if(image .ne. 1 .or. image .ne. nimFPI(1)) noldFPIat=0
cc          do ia=1,natomsold
cc            if(itag(ia) .ge. 100) then
cc                if(image .eq. 1 .or. image .eq. nimFPI(1)) then
cc                  itag(ia)=1
cc                else
cc                  itag(ia)=0
cc                  noldFPIat=noldFPIat+1
cc                  ioldFPIat(noldFPIat)=ia
cc                endif
cc           endif
cc          enddo
cc         WRITE(lunout,288) noldFPIat,nfatms

          nFPI=0

c           the total number of atoms should now be nFPI*(nimFPI(1)-1) less               
c           than before, check this:                                                      
          if(itot .ne. natoms) then
            write(lunout,*) ' ** ERROR loop:   itot .ne. natoms '
            write(lunout,*)'itot=',itot,'nFPI=',nFPI,'nimFPI=',nimFPI(1)
            write(lunout,*)'    natoms = ',natoms,'  natomsold=',natomsold
          endif

          indupd=1
c          make sure a list of interactions is created now.                               
            
	  numIMAGE=1

          CALL FORCE
c$$$          if (genConstraints .and. iqkmin) call projectForce
          if (genConstraints) call projectForce
c          CALL KINET

c       Find the RMS and maximum force on the atoms in this image:
          fsum2=0.0
          fmax=0.0
          nfatms=0
          do iafc=1,noldFPIat
c          should only include movable atoms that used to be in chains:
              iaf=ioldFPIat(iafc)
              ind=3*(iaf-1)
              fsum2=fsum2+fa(ind+1)**2+fa(ind+2)**2+fa(ind+3)**2
              if(abs(fa(ind+1)) .gt. fmax) fmax=abs(fa(ind+1))
              if(abs(fa(ind+2)) .gt. fmax) fmax=abs(fa(ind+2))
              if(abs(fa(ind+3)) .gt. fmax) fmax=abs(fa(ind+3))
          enddo

          RMSforce=fsum2/(3.0*noldFPIat)

          if(image .eq. 1) refimutot=utot
c              If few atoms are included as chains, the second image can be lower
c              in energy than the first, and then it is better to use the second
c              image as baseline:
          if(image .eq. 2 .and. utot .lt. refimutot) then
             sliding=utot-refimutot
             refimutot=utot
             difpot=utot-refimutot
             write(lunout,603) image,UTOT,difpot,ttot,RMSforce,fmax
             write(lunout,605) sliding
605          format('#f  ',23x,' sliding down by ',g12.4,
     +              ',  shift ref to 2nd im.')
          else
            difpot=utot-refimutot
            write(lunout,603) image,UTOT,difpot,ttot,RMSforce,fmax
603         format('#f  ',i5,g25.9,g12.4,g15.5,2f10.6)
          endif
          relenimage(image)=difpot

c      Dump the coordinates of this image into file coFPI.con:
cc          call conout(luncoFPI,-1,ax)
c      Dump coordinates of the maximum energy image into file cmaxim.con
            if(Inudge .and. (utot .gt. enmax)) then
              enmax=utot
              nummaxim=image
              open(lunmaxim,file='cmaxim.con')
              call conout(lunmaxim,-1,ax)
              close (lunmaxim)
            endif

c      Write the displacement from previous image to this image for 
c            the first 5 FPI chains:
c            write into file outFPI.dat

cc         ntimesd=0.2*nFPIold+1
         ntimesd=1
         ifirstch=-4
         ilastch=0
         do ich=1,ntimesd
           ifirstch=ifirstch+5
           ilastch=ilastch+5
           if(ilastch .gt. nFPIold) ilastch=nFPIold
cERB           write(lunoutFPI,604) 
cERB     +            image,(ax*displim(ipp),ipp=ifirstch,ilastch),difpot
604        format('#d  ',i5,10g12.4,g15.6)
         enddo

600       continue
c              close loop over images.

cH2O------------
       else         ! if not all in chain option do above, else below

          dl = 0.
          do image = 1, nimFPI(1)
             utotim(image) = utotim(image) * nimFPI(1)
             if (image .eq. 1) then
                dl2 = 0.d0
             else
                dl2 = 0.
                do i = 1, (nHydrogens+nOxygens)/nimFPI(1)
                   indAt = 3 * nimFPI(1) * (i-1) + 3 * (image - 1)
                   do j = 1, 3
                      ind = indAt + j
                      dl1 = ra(ind)-ra(ind-3)
                      call pbc(dl1, box(j)/ax)
                      dl2 = dl2 + dl1*dl1
                   end do
                end do
             end if
             dl = dl + sqrt(dl2) * ax


             if(image .eq. 1) refimutot=utotim(image)
             if(image .eq. 2 .and. utot .lt. refimutot) then
                sliding=utotim(image)-refimutot
                refimutot=utotim(image)
                difpot=utot-refimutot
cERB                write(lunout,603) image,UTOTim(image),difpot,ttot,
cERB     $               RMSforce,fmax
                write(lunout,663) image, dl, UTOTim(image), difpot,
     $               t(image), sqrt(fuperpT2(image))*ax, 
     $               sqrt(faT2(image))*ax
                write(lunout,605) sliding
             else
                difpot=utotim(image)-refimutot
cERB                write(lunout,603) image,UTOTim(image),difpot,ttot,
cERB     $               RMSforce,fmax
                write(lunout,663) image, dl, UTOTim(image), difpot,
     $               t(image), sqrt(fuperpT2(image))*ax, 
     $               sqrt(faT2(image))*ax
             endif
          end do

       end if

663         format('#f  ',i4, f10.6, g20.9,g12.4,g15.5,2f10.6)


cAllinChain
cERB          end if
cAllinChain

	numIMAGE=nimFPI(1)
cc      endif

c   Get optimal estimate for barrier height if this is an NEB calc:
      if(Inudge) then
c       get error bounds for energy barrier based on taking highest
c       energy image as minimum estimate and increment that by
c       the difference between highest and next highest energy image
c       to get the maximum estimate:
          barrmin=relenimage(nummaxim)
c       get next highest energy
          ennextmax=relenimage(nummaxim-1)
          if(relenimage(nummaxim+1) .gt. ennextmax) 
     +              ennextmax=relenimage(nummaxim+1)
          halfincrement=0.5*(relenimage(nummaxim)-ennextmax)
          barrave=barrmin+halfincrement
          barrmax=barrmin+2.0*halfincrement

         WRITE(lunout,289) barrave,halfincrement
289      format('#'/'# NEB:  Barrier height = ',g12.4,'  +/- ',f7.3)

      endif
c        close  `if(Inudge) ...'

         WRITE(lunout,288) noldFPIat
288      format(/'#  noldFPIat: # of tag 0 atoms from chains =',i7)

      endif
c                    close `if(nFPI .gt. 0) then'                                         

c   ----------------------------------------------------------------------                
c                                                                                         
c  For CPU timing:                                                                        
cG77      CALL SECOND(cputal)
      cputil=cputal-cputbl
c                                                                                         
c                                                                                         
      WRITE(lunout,280) maxneb
280   format('#'/'# Maximum number of interacting neighbors = ',i10)
c                                                                                         
      write(lunout,282) maxnpr
282   format('# Maximum number of pairs =           ',i15)
c                                                                                         
      WRITE(lunout,281) nstupd
281   format('#'/'# At the end, nstupd = ',i6)
c                                                                                         
      WRITE(lunout,283) nstochcoll,nmasscoll
283   format('#'/'# Total number of stoch coll on individual atoms, ',
     +                'nstochcoll= ',i6/
     + '#   -     -     - massive collisions,        nmasscoll = ',i6)

c        --------------------------------------------------------------------             

cG77      CALL SECOND(cputf)
cG77      CPUTIM=(cputf-cputbl)

c666      call etime(cputf2)
cG77      WRITE(lunout,1295) CPUTIM
      WRITE(lunout,1295) cputf2(1)-cputbl2(1), cputf2(2)-cputbl2(2)
cG77 1295  FORMAT('#'/'#'/'#  TOTAL CPU TIME in loop:    ',g12.4,' seconds')

      print *, ' ********  FIN *********'

      do i = 1, 8
         do j = 1, 3
            r1(j) = (ra(3*(2*i-2)+j) + ra(3*(2*i-1)+j) + 16.d0 * 
     $           ra(3*(16+i-1)+j)) / 18.d0 * ax
         end do
         write(98,'(3f15.5)') r1(1), r1(2), r1(3)
      end do


      if (1 .lt. 0) then
         CALL FORCE

         open(88, file='coord',status='unknown')
         open(98, file='force',status='unknown')

         write(88, '(A)') '=========================================='
         do i = 1, 288
            write(88, '(3f16.9,2i6)') ra(3*(i-1)+1)*ax, ra(3*(i-1)+2)*ax
     $           ,ra(3*(i-1)+3)*ax, 0, i
         end do
         write(88, '(A)') '=========================================='
         do i = 1, 288
            write(88, '(3f16.9,2i6)') va(3*(i-1)+1)*ax, va(3*(i-1)+2)*ax
     $           ,va(3*(i-1)+3)*ax, 0, i
         end do
         
         write(98, '(A)') '=========================================='
         do i = 1, 288
            write(98, '(3f16.9,2i6)') fa(3*(i-1)+1), fa(3*(i-1)+2), fa(3
     $           *(i-1)+3)
         end do
      end if
      close(88)
      close(98)

1295  FORMAT('#'/'#'/'#  TOTAL CPU TIME in loop: '/
     $     '#                     user time ', g10.4, ' seconds'/
     $     '#                   system time ', g10.4, ' seconds'/'#')
c                                                                                         
cc      if(cputim .lt. 0.1e-10) go to 999                                                 
cc      WRITE(lunout,1297) CPUTIG,cputig/cputim*100.                                      
cc1297   format('             cpu time in gagafe:   ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1291) CPUDEF,cpudef/cputim*100.                                      
cc1291   format('                  300  calc. DEL   ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1292) CPUTPB,cputpb/cputim*100.                                      
cc1292   format('                  320  apply PBC   ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1293) CPUPBS,cpupbs/cputim*100.                                      
cc1293   format('                       skewed PBC  ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1299) CPUTR2,cputr2/cputim*100.                                      
cc1299   format('                  371  calc. R2    ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1290) CPU372,cpu372/cputim*100.                                      
cc1290   format('                  372  make list   ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1284) CPU380,cpu380/cputim*100.                                      
cc1284   format('                  380  NoUpd delpr ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1285) CPU381,cpu381/cputim*100.                                      
cc1285   format('                  381  NoUpd PBC   ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1286) CPU382,cpu382/cputim*100.                                      
cc1286   format('                  382  NoUpd R2    ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1583) CPU383,cpu383/cputim*100.                                      
cc1583   format('                  383  Fix large R2',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1570) CPU370,cpu370/cputim*100.                                      
cc1570   format('                  370  calc. force ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1577) CPU377,cpu377/cputim*100.                                      
cc1577   format('                  377  zero force  ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1576) CPU376,cpu376/cputim*100.                                      
cc1576   format('                  376  copy force1 ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1578) CPU3761,cpu3761/cputim*100.                                    
cc1578   format('                  3761  virial     ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
cc      WRITE(lunout,1575) CPU375,cpu375/cputim*100.                                      
cc1575   format('                  375  copy force2 ',g12.4,'   seconds',                 
cc     +     f9.3,'%')                                                                    
c                                                                                         
c999   continue

      RETURN
      END

c-----------------------------------------------------------------------
 
      subroutine pbc(r, a)

      real*8 r, a

      if(r .gt. a/2.0d0) then
         r = r - a
      elseif(r .lt. -a/2.0d0) then
         r = r + a
      end if

      return
      end
