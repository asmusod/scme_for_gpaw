c   EAM   version c:    make tags 100-199 indicate FPI with fixed ends.                   
c                       and spring constants depend on type only (not for each            
c                       FPI chain separately).                                            
c                     Use common blocks in ../commonblks/..                               

c       3dmdD:     version D, (in development) needs one more line in input that          
c                           has option for dumping data on viscosity, thermal cond.       
c       3dmdC:      version C,  needs one more line in input that                         
c                              specifies whether any parameters are changing.             
c       3dmdB:      version B,  Feb 13.                                                   
c       3dmdA:      version A, Jan20 1990.                                                
c                                                                                         
c     MAIN PROGRAM FOR MOLECULAR DYNAMICS OF 3-D SYSTEMS.                                 
c         Canonical, isobaric, microcanonical, isobaric-isoenthalpic.                     
c       modify version A of SURF, oct 89.                                                 
c                                                                                         
      subroutine main(natm, aseCoords,aseBox,eQM,FAOUT,EPOTOUT,ETOUT)
      implicit real*8 (a-h,o-z)
      integer natm
      real*8 aseCoords(3*natm)
      dimension FAOUT(3*natm)
      real*8 ETOUT(3,natm / 3)
      real*8 EPOTOUT
      real*8 aseBox(3)
      real*8 eQM(3,natm/3)

c     since this is stricly water only pot. #Atomtypes =2 always
c      logical Iqkmin, IPRES, ILINCH Inudge
cf2py integer intent(in) :: natm
cf2py real(8) intent(in) :: aseCoords 
cf2py real(8) intent(in) :: aseBox
cf2py real(8) intent(in) :: eQM
cf2py real(8) intent(out) :: FAOUT
cf2py real(8) intent(out) :: EPOTOUT 
cf2py real(8) intent(out) :: ETOUT 

c
       include '../commonblks/parameters.cmn'
       include '../commonblks/comconf.cmn'
       include '../commonblks/comtime.cmn'
c      combaths has irigidmolecules
       include '../commonblks/combaths.cmn'
       include '../commonblks/comluns.cmn'
       include '../commonblks/comgeom.cmn'
       include '../commonblks/constraints.cmn'

       print *, eQM

c  Hardcode box dims to make stuff work before 
c  inputting from python
c  box dimensions 
      ax = aseBox(1)
      ay = aseBox(2)
      az = aseBox(3)

      alpha=90.
      beta=90.
      gamma=90.
c  Thermostat params (needs to be there but not used)
      Pmass = 1.
      Volm = 0.0
c  Atom type masses
      AMASS(1) = 1.
      AMASS(2) = 16.
c  Number of atom types - water only: always 2
      NATYPE = 2

c  Stuff from READIN - we don't want no inp.dat
c    MD times: keep hardcoded, so we only
c              initialize and thus get forces
      tottim = 1.
      stpsz  = 1. 
c  more hardcoded stuff from readin
c  proobably not necessary - check at some point
      Iqkmin = .FALSE.
      tolforce = 0.1e-10
      IPRES = .FALSE.
      ILINCH = .FALSE.
      Inudge = .FALSE.
      IrigidMolecules = .TRUE.

c                                                                                         
c   variables for logical unit numbers:                                                   
      lunout=10
      luninp=16
      lunci=17
      lunco=11
      luncm=12
      lunins=13
      luntgr=14
      lunaco=18
      lunvis=19
      lunattr=20
      luncoFPI=21
      lunoutFPI=22
      lunqckmin=23
      lunoutattr=24
      lunmaxim=25

CAA Call potinit to get some parameters necessary for the water potential
      CALL potinit()
CAA

c  There are two input files:                                                             
      open(luninp,file='inp.dat')
      open(lunci,file='ci.con')
c   and there can be six output files: (only open one here, the others are                
c                                       opened in subroutine LOOP)                        
      open(lunout,file='scme_out.dat')
c      write(6,*) '  Get the coordinates from ASE:'
      NATOMS = natm
c     1 = HYDROGEN FFS!
c     2 = OXYGEN
      NATMS(1) = NATOMS * 2 / 3
      NATMS(2) = NATOMS / 3
      do i = 1,3*natm
         RA(i) = aseCoords(i)
c         print *, RA(i)
      end do

c      CALL CONIN
      close(lunci)
c        reads the configuration file 'ci'.                                               
c                                                                                         
cc      write(lunout,203)                                                                 
cc203   format(/'   after calling conin ...')                                             
c      write(6,*) '  Read in the inp.dat file - NO'
c      CALL READIN
      close(luninp)
c        reads the input file on unit 7 to get the job control parameters,                
c        writes some of the control parameters to unit 5.                                 
c                                                                                         
cc      write(lunout,200)                                                                 
cc200   format(/'   after calling readin ...')                                            
c                                                                                         
c      write(6,*) '  Call Setup:'
      CALL SETUP(eQM,FAOUT,EPOTOUT,ETOUT,natm)
c        prepares the input configuration for dynamics, initializes etc.                  
c                                                                                         
cc      write(lunout,205)                                                                 
cc205   format(/'   after calling setup..')                                               
c      write(6,*) '  Call Loop:'                                                                                         
c      CALL LOOP
c        performs the dynamics for the required number of timesteps.                      
c                                                                                         
      end subroutine main
      END


