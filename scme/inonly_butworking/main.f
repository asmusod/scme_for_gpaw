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
      subroutine main(natms, FAOUT)
!f2py integer intent(in) :: natms
!f2py real(8) intent(out) :: FAOUT

      implicit real*8 (a-h,o-z)
c
      
       include '../commonblks/comluns.cmn'
      dimension FAOUT(3*natms)
      intent(out) FAOUT
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
      open(lunout,file='out.dat')
      write(6,*) '  Read in the coordinates:'
      CALL CONIN
      close(lunci)
c        reads the configuration file 'ci'.                                               
c                                                                                         
cc      write(lunout,203)                                                                 
cc203   format(/'   after calling conin ...')                                             
      write(6,*) '  Read in the inp.dat file:'
      CALL READIN
      close(luninp)
c        reads the input file on unit 7 to get the job control parameters,                
c        writes some of the control parameters to unit 5.                                 
c                                                                                         
cc      write(lunout,200)                                                                 
cc200   format(/'   after calling readin ...')                                            
c                                                                                         
      write(6,*) '  Call Setup:'
      CALL SETUP(FAOUT)
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


