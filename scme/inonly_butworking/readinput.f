c   EAM version e:  Add option for Quickmin by selecting negative collision rate.
c   EAM version c:  Spring const for each type only (not each FPI chain)                  
c   Feb 25. 92:  add indf1 in input and calculate iorderf                                 
c   EAM version b91: add FPI.  If FPI atoms in ci.con then read in line with              
c                              spring constant.                                           
C        July 3. 1990:     add viscosity output option.                                   

C     THIS SUBROUTINE READS THE INPUT FILE:                                               
C                                                                                         
      SUBROUTINE READIN
c                                                                                         
      implicit real*8 (a-h,o-z)

      character header(80)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comoutcntrl.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comtgr.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/compotent.cmn'
      include '../commonblks/comdeposit.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comintlis.cmn'
      include '../commonblks/constraints.cmn'

      common /cputim/ cput0,cputas,cputal,cputig,cpudel,cputr2,cputbl
     +               ,cpudef,cputpb,cpupbs,cpucpr,cpu372,cpu370,cpu373
     +               ,cpu380,cpu381,cpu382,cpu383,cpu375,cpu377,cpu376
c         these variables are used for finding CPU time in various parts.                 

      data pi /3.14159265358979312d0/
      data degperrad / 57.2957795130823229d0 /

C CALL TIMER TO ANALYSE CPU TIME:                                                         
c (The subroutine SECOND is defined here in this code and is macine dependent).           
cG77      CALL SECOND(cput0)
      cputig=0.
      cpudel=0.
      cputr2=0.
      cputbl=0.
      cpudef=0.
      cputpb=0.
      cpupbs=0.
      cpu372=0.
      cpu370=0.
      cpu373=0.
      cpu375=0.
      cpu376=0.
      cpu377=0.
      cpu380=0.
      cpu381=0.
      cpu382=0.
      cpu383=0.
C                                                                                         

c   -----------------------------------------------------------------------

      read(luninp,10) (header(i),i=1,80)

C  WRITE TO OUTPUT FILE:                                                                  

      WRITE(lunout,190)
190   FORMAT('#'/'#'/
     +    '# MOLECULAR DYNAMICS in 3-D, md3dEAMf code  June 1994'/'#')

      write(lunout,11) (header(ii),ii=1,78)
11    format('#',78a1)

      if(az .le. 0.99e-06) write(lunout,191) az
191   format('#'/
     +       '# This is a two dimensional simulation since az = ',g12.4)

c           --------------------------------------------------------------

C   READ PROCESS INFORMATION:                                                             

      READ(luninp,10)
c                skip second line.                                                        
10    FORMAT(80A1)
c                                                                                         
      read(luninp,*) tottim, stpsz
c                total time of simulation.                                                
c                STPSZ is the time steplength,                                             
      if(stpsz .le. 0.00) go to 999
      ntotstps=tottim/stpsz+0.5

           WRITE(lunout,300) STPSZ,TOTTIM,ntotstps
300        FORMAT('# TIME STEP SIZE        = ',f12.4/
     +          '# TOTAL SIMULATION TIME = ',f12.4/
     +          '# NUMBER OF TIME STEPS  = ',i12)                                                                                         
                                                                                        

      read(luninp,*) incol,temmas
c                 Do initial massive collision?                                           
c                 temperature for  - || -     .    

           IF(INCOL) WRITE(lunout,210) temmas
210        FORMAT('#'/'#   DO INITIAL MASSIVE COLLISION:'/
     +           '#      with Temp = ',g12.4)
c                                                                                        
                                       
      read(luninp,*) ITEMP,imasv,textl,colrat
c                 itemp:   keep temperature constant?                                     
c                 imasv:   massive collisions?  (otherwise single atom per                
c                                                collision)                               
c                 textl:   temperature                                                    
c                 colrat:  number of collisions per unit time (per TAU)                   
        crate=colrat*stpsz
c               'crate' is the number of collisions per timestep.                         
        if(itemp) then
         if(imasv) then
c         Do periodic massive collisions:                                                 
            irate=.false.
            nmasv=1./crate+0.5
c                       number of timesteps between massive collisions.                   
            if(nmasv .le. 0) then
               nmasv=1
               write(lunout,*)  
     +         '#  **  WARNING readin: Cannot have more massive',
     +         '  collisions than timesteps.'
               write(lunout,*)'#    Change to collision every timestep.'
            endif
            crate=0.0
         else
c         Do individual stochastic collisions:                                            
            irate=.true.
            nmasv=9999999
         endif
      endif

c   write:
          IF(ITEMP) then
            IF(imasv) WRITE(lunout,220) colrat
220         FORMAT('#'/'#  CONSTANT TEMPERATURE DYNAMICS:'/
     +              '#    USING MASSIVE COLLISIONS,',f9.3,' per TAU')
            IF(.not. imasv) WRITE(lunout,221) colrat
221         FORMAT('#   CONSTANT TEMPERATURE DYNAMICS:'/
     +          '#      USING INDIVIDUAL STOCH. COLLISIONS,',
     +                  f9.3,' per TAU')
            WRITE(lunout,301) TEXTL
301         FORMAT('#     BATH TEMPERATURE     =',F10.5)
          endif

c                              -------

      read(luninp,*) Iqkmin, tolforce
c         Check if  quick minimization is to be used, and tolerance for stopping:
c         Quick minimization zeroes velocity perpendicular to force at each time
c         step, and if the parallel velocity is in opposite direction as
c         force, it is zeroed too.  
c         tolforce:    the run will be stopped when the RMS value of the force
c                      has dropped below `tolforce' if that happens within
c                      the total time of the run which is given in the third
c                      line of inp.dat.
         if(Iqkmin .and. Ipres) then
           write(lunout,310)
310        format(/' *** You should not select quick minimization when',
     +             ' doing constant pressure, to be added later. ')
           stop
         endif 

         if(Iqkmin) write(lunout,307) 
307       format(/'#  Quickmin, zero velocity perpendicular to force')

         if(Iqkmin .and. Itemp) then
          write(lunout,311) 
          write(6,311) 
311       format(/'#  *** Do not couple to heat bath while doing ',
     +     'quick minimization.'/'#   Heat bath will be turned off.'/)
          Itemp=.false.
         endif

c                              -------

      IPCOL=.false.
      read(luninp,*) IPRES,pextl,plambd
c                 IPRES:   keep pressure constant?                                        
c                 pextl:   external pressure                                              
c                 plambd:  weight of stochastic component in piston velocity.             
c                                                                                         
      if(abs(plambd) .gt. 0.1e-10 .and. IPRES) IPCOL=.true.
c                 IPCOL is .true. when piston gets hit by stochastic collisions           

          IF(IPRES) then
              WRITE(lunout,222) pextl,plambd
222           FORMAT('#'/'#   CONSTANT PRESSURE DYNAMICS:'/
     +            '#           external pressure    = ',g12.4/
     +            '#           weight of stoch coll = ',g12.4)
          endif

      read(luninp,*) ILINCH,indpara,ratlinch
c                 ILINCH:   change some parameter linearly with time?                     
c                 indpara:  indicates which parameter should be changed                   
c                        = 1  for temperature                                             
c                        = 2  for pressure                                                
c                 ratlinch:  rate of change per tau.                                      

         IF(ILINCH) then
             if(indpara .eq. 1) WRITE(lunout,223) ratlinch
223          FORMAT('#'/'#   CHANGE BATH TEMPERATURE LINEARLY:'/
     +            '#           rate (change per tau)    = ',g12.4)
             if(indpara .eq. 2) WRITE(lunout,224) ratlinch
224          FORMAT('#'/'#   CHANGE EXTERNAL PRESSURE LINEARLY:'/
     +            '#           rate (change per tau)    = ',g12.4)
         endif

c                    ------------------------------------                                 

cFPI:     Read spring constants for FPI chains:                                           
      read(luninp,*) (sprconFPI(k),k=1,natype)

      read(luninp,*) Inudge, phiswitchsmall, phiswitchlarge
c         Check if nudging should be used and calculate 
c         cos(angles)
         if(Inudge) then
            cosPhismall=cos(phiswitchsmall/degperrad)
            cosPhilarge=cos(phiswitchlarge/degperrad)
c          Nudging throws out certain components of the force,
c            then SADDLEFIND is called instead of FPIforce.
            write(lunout,306) phiswitchsmall,phiswitchlarge
306         format(/'#  Find Minimum Energy Path using nudged elastic ',
     +      'band method'/'#  zero some force components of the chain.'/
     +      '#  use switching function when phi is between',2g10.2)
            if(phiswitchsmall .gt. phiswitchlarge) then
             write(lunout,312)
312          format('  Cannot have phiswitchsmall >',
     +              ' phiswitchlarge')
             stop
            endif
         endif
 
c    ------------------------------------------------------------------------------

cDeposit:                    
c  Read parameters for Atom Deposition:
      read(luninp,*) Ideposit,dz0vd,v0vd,taftvd,tbtwmc,
     +               tbfnex,afract,ztcontrl
c               if Ideposit -  if true then atoms are deposited (ATOM DEPOSITION),
c                  dz0vd  -  is the initial z coordinate of the VD
c                           atoms, above the highest interacting surface atom
c                           (zsurf)
c                  v0vd  -  is the initial z velocity of the atoms,
c                  taftvd - is the time for dynamics after starting the VD atom
c                           and until a massive stochastic coll. is performed.
c                  tbtwmc - is the time for dynamics between the two massive 
c                           stochastic collisions. 
c                  tbfnex - is the simulation time after the second massive
c                           collision and before the next atom deposition. 
c                  afract - is the fraction of A atoms in the VD atoms.
c                  ztcontrl - maximum z value for atoms that get hit with
c                             stochastic collisions.
c
      if(Ideposit) then
c           calculate the number of steps between the various AD events.
        naftvd=taftvd/stpsz+0.5
c           the number of timesteps after an atom is sent towards the surface
c           until the first massive stochastic collision is performed.
        nbtwmc=tbtwmc/stpsz+0.5
c           number of time steps between the two massive collisions
        nbfnex=tbfnex/stpsz+0.5
	if(nbfnex .le. 0) nbfnex=1
c           number of time steps after the second massive
c           collision until the next atom is sent towards the surface,
c           make sure that the massive collision is
c           done at least one time step before new AD atom is sent in.
        ntotvd=naftvd+nbtwmc+nbfnex
c           the total number of steps between depositing atoms.
        ttotvd=ntotvd*stpsz
c           total simulation time per each deposited atom.
      endif

          if(Ideposit) write(lunout,229) dz0vd,v0vd,taftvd,naftvd,
     +          tbtwmc,nbtwmc,tbfnex,nbfnex,ttotvd,ntotvd,afract,ztcontrl
229       format('#'/'#  ATOM DEPOSITION:'/
     +        '#    atoms start ',g12.4,' above zsurf'/
     +        '#    with velocity Vz0 = ',g12.4/
     +    '#    time before 1st mass. coll =',g12.4,' # of steps =',i5/
     +    '#    time between massive coll. =',g12.4,' # of steps =',i5/
     +    '#    time after 2nd mass. coll. =',g12.4,' # of steps =',i5/
     +    '#    total time per atom depos. =',g12.4,' # of steps =',i5/
     +    '#    fraction of A atoms in AD  =',g12.4/
     +    '     maximum z value for atoms hit by stoch coll. = ',g12.4)


c        ------------------------------------------------------
c-----------------------------------------------------------------H2O /*
c     If IrigidMolecules=.true. the rattle algorithm will be used
c     assuming that the molecules are H2O. ERB 29/9/94
      read(luninp,*) IrigidMolecules

      print '(A, $)', 'Are molecules rigid? '
      print *, IrigidMolecules
      if (IrigidMolecules) then
         print *
         print '(A)', ' Water molecules will be considered rigid bodies'
         print '(A)', ' with the parameters input in inp.dat file.'
         print *
      end if

c     Are we imposing generalized constraints of the form A*v=0?
c     If true, the program will look for the vector A in file 
c     constraints.dat. The file has the following form:
c     A1        A2        A3         1
c     A4        A5        A6         2
c
c     ......................
c
c     etc.
c     There is no label for the 2nd component!!

      read(luninp,*) genConstraints
c-----------------------------------------------------------------H2O */
               
c-----------------------------------------------------------------------

                                                               
C  READ CONTROL PARAMETERS FOR OUTPUT:                                                    
                                                                                        
      READ(luninp,10)
      READ(luninp,10)
      READ(luninp,*) ngrp
c                     NGRP is the # of steps to inlcude in averages.                      

      write(lunout,304) ngrp
304   format(/'#   Number of steps to inlcude in averages, ngrp =',i5)

cw      write(6,*) '   readin:  ngrp = ',ngrp

      READ(luninp,*) IDUMP,NDUMP,NWAIT,indvel,zminwr
c                     if IDUMP = .true. then intermediate configurations are
c                     dumped every NDUMP steps starting after NWAIT steps. 
c                     if indvel=-1 then velocities are not written out.   
c                     zminwr = minimum z value for including atoms in output

C                                                                       
      READ(luninp,*) ITEMPGR,ngrlx0,ngrly0,ngrlz0
c             'ngrlx'  is the number of meshes in the x direction, etc.
      if(ITEMPGR) then
         ngrlx=(ngrlx0+1)/2+(ngrlx0-1)/2
         ngrly=(ngrly0+1)/2+(ngrly0-1)/2
         ngrlz=(ngrlz0+1)/2+(ngrlz0-1)/2
         if(ngrlx .ne. ngrlx0) write(lunout,*) 
     +              '**  WARNING: ngrlx must beodd'
         if(ngrly .ne. ngrly0) write(lunout,*) 
     +              '**  WARNING: ngrly must beodd'
         if(ngrlz .ne. ngrlz0) write(lunout,*) 
     +              '**  WARNING: ngrlz must beodd'
      endif
      IGRP=.true.
      IRMS=.true.

          IF(ITEMPGR) write(lunout,302) ngrlx,ngrly,ngrlz
302       format('#'/
     +           '#       Temperature grid:     mesh size x:   ',i5/
     +           '#                                       y:   ',i5/
     +           '#                                       z:   ',i5)

                                                                                         
      READ(luninp,*) INSTA,NINST
c                     if IINST = .true. then instantaneous values are written             
c                     to unit lunins.                                                     
        if(INSTA) write(lunout,305) NINST
305     format('#  Number of timesteps betw. instantaneous values =',i5)

      READ(luninp,*)  IAVECOO,nstaveRA
c                nstaveRA is the number of time steps to use in averaging                 
c                         the coordinates of the atoms.                                   

        if(IAVECOO) write(lunout,303) nstaveRA
303     format('#  Number of timesteps in coord. aver., NstaveRA =',i5)
          
                                                                               
      READ(luninp,*) IVISC,NVISC,NWVISC
c                if IVISC = .true. then data                                         
c                that can be used to calculate viscosity and thermal conductivity    
c                are dumped every NVISC steps starting after NWVISC steps.           

        if(IVISC) write(lunout,308) NVISC,NWVISC
308     format('#  Write info on viscosity, NVISC,NWVISC =',2i7)
 

      READ(luninp,*) Iattraj,indatwr,nstatwr
c                if Iattraj = .true. then coordinates, force and energy of atom
c                   indatwr  is written out to file `attraj.dat' every
c                   nstatwr steps.

        if(Iattraj) write(lunout,309) indatwr,nstatwr
309     format('#   Write trajectory of atom',i4,' every ',i6,' steps')
 
c  ----------------------------------------------------------------------------                                                                                         
c  ----------------------------------------------------------------------------                                                                                         
c                                                                                         
C  READ PARAMETERS FOR THE POTENTIAL:  (only the gagafe routine knows how to              
c                  interpret these parameters. They are simply passed on to the           
c                  common block 'POTENT' by this routine.)                                
C                                                                                         
      READ(luninp,10)
c                         skip line.
      READ(luninp,10)
c                         header for Potential Parameters.
                                                                                      
c   The number of different interactions (A-A, B-B, A-B, etc.) is                         
      ninteract=NATYPE*(NATYPE+1)/2
c         where NATYPE is the number of components.                                       
CAA
CAAcw      write(6,*) '  readin:  ninteract = ',ninteract
CAA
CAA      read(luninp,*) (RCUT(ii),ii=1,ninteract)
CAA
CAAcw      write(6,*) '  readin: RCUT(1-3) = ',RCUT(1),RCUT(2),RCUT(3)
CAA
CAA      read(luninp,*) (Rskin(ii),ii=1,ninteract)
CAAc                        rcut is the range of the potential.                              
CAAc                        rskin is the range of the buffer range.                          
CAAc                        The order for these parameters is: AA, BB, AB.                   
CAA
CAAcw      write(6,*) '  readin:  rskin(1-3)= ',(rskin(i),i=1,3)
CAA
CAA      read(luninp,*) indf1
CAAc              indf1 points to the parameter that corresponds to first order in           
CAAc                    rho in the polynomial in f(rho).                                     
CAA      npotpar=0
CAA      do ip=1,MAXPOTPAR
CAA         read(luninp,*,end=99,err=98) (potpar(ip,ii),ii=1,ninteract)
CAAcc     +                    ,(potheader(k,ip),k=1,30)                                     
CAAcc         write(6,*) '   READIN:  read parameter # ip = ',ip                             
CAA         npotpar=npotpar+1
CAA      enddo
CAAcc361   format(3g20.0,30a1)                                                               
CAA
CAA      write(lunout,*) '  ERROR READIN:  only ',MAXPOTPAR,
CAA     +               ' potntl. parameters  are allowed'
CAA      stop
CAA99    continue
CAA      write(6,*) '   READIN:  end  in read,   # of pot par = ',(ip-1)
CAA      go to 97
CAA98    continue
CAA      write(6,*) '   READIN:  err  in read,   # of pot par = ',(ip-1)
CAA97    continue

c   calculate the order of the polynomial for the embedding function:                     
      iorderf=npotpar-indf1+1

CAA      do 259 i=1,3
CAA         rcut2(i)=rcut(i)**2
CAA         rskin2(i)=rskin(i)**2
CAA259      continue
c                                                                                         
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                 
c     ---------------------------------------------------------------------


cc      WRITE(lunout,362)
362   format('#  ---------------------------------------------')



      RETURN
999   continue
      write(lunout,9999) stpsz
9999  format('#'/'**  ERROR 999 in subr. READIN,   stpsz = ',g12.4)
      stop
      END



