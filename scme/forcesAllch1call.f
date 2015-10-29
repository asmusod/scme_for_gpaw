c   ForcesAllChains.f     for use when all atoms are represented as chains
c                         then this routine selects out coordinates corresponding
c                         to each image and sends to GAGAFE.
c   version f:
c   version e93:   add gregs SPF routine for chain relaxation
c   EAM version b91:   add FPI forces and potential by calling FPIforce                   
c                                                                                         
C     THIS SUBROUTINE COMPUTES THE FORCES by calling GAGAFE for each type of              
c     interaction.                                                                        

      SUBROUTINE FORCE(eQM,FAOUT,ETOUT,natm)

      implicit real*8 (a-h,o-z)
      integer natm
      real*8 eQM(3,natm/3)
      intent(in) eQM
      dimension FAOUT(*), ETOUT(*)
      intent(out) FAOUT, ETOUT

      include '../commonblks/parameters.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comenergy.cmn'
      include '../commonblks/comenperat.cmn'
      include '../commonblks/comluns.cmn'

      common /unscale/RALOCAL(MAXCOO)

cc      equivalence (Potenpat(1),Potpat(1))                                               

      dimension vir(6)
      dimension natmsold(MAXCOMP),FALOCAL(MAXCOO)
c$$$      dimension utotim(MAXATOMS)

      logical  indH2O
      integer itagl(MAXATOMS)

c     -----------------------------------------------------------------------

      indH2O=.true.

      DO I=1,3*NATOMS
        FA(I)=0.
        FAOUT(I) = 0.
        FAlocal(i)=0.0
        RALOCAL(I) = 0.0
c             the array FA will contain the force on each atom after all calls            
c             to GAGAFE have been completed.                                              
      enddo
      do i=1,natoms
        Potenpat(i)=0.
c           the array Potenpat will contain the potential per atom after all calls        
c           to GAGAFE have been completed.                                                
      enddo
c   -------------------------------------------------------------------


c   Extract successive images and send each one to GAGAFE:

      if(nFPI .gt. 0) then
c      make sure ALL atoms are in fact in FPI chains:
         if(natoms .ne. nFPI*nimFPI(1)) then
           write(lunout,*)'ERROR, forcesAllChains: natoms=',natoms
           write(lunout,*)'    while nFPI*nimFPI(1)=',nFPI*nimFPI(1)
           stop
         endif

        natomsold=natoms
        natoms=nFPI
        do it=1,natype
          natmsold(it)=natms(it)
          natms(it)=nFPIcom(it)
        enddo
      endif
c                    close `if(nFPI .gt. 0) then'

      if(nFPI .eq. 0) nimFPI(1)=1
c              fake this temporarily to get the case of no FPIs to run.

      if (nimFPI(1) .eq. 1) then
         do na=1,natoms
            itagl(na) = itag( na )
         end do
      else
         do na=1,natoms
            itagl(na) = itag( nimFPI(1) * (na-1) + 2 )
         end do
      end if

      do 250 im=1,nimFPI(1)

cw          write(*,*) ' '
cw          write(*,*) ' From forcesAllChains:   time= ',time
cw          write(*,*) '   coords of image ',im

          do na=1,natoms
             ind1=3*(na-1)
             ind2=3*((na-1)*nimFPI(1)+im-1)
             ralocal(ind1+1)=ra(ind2+1)*ax
             ralocal(ind1+2)=ra(ind2+2)*ax
             ralocal(ind1+3)=ra(ind2+3)*ax


cw             x=ralocal(ind1+1)
cw             y=ralocal(ind1+2)
cw             z=ralocal(ind1+3)
cw             write(*,*)' Ra=',x,y,z

             FAlocal(ind1+1) = 0.0
             FAlocal(ind1+2) = 0.0
             FAlocal(ind1+3) = 0.0

          enddo

cw        write(lunout,*) '   Force:     Enter the routine,  NATYPE=',NATYPE              
cw        write(lunout,*)RA(1),RA(2),RA(3),ax                                             
cw        write(lunout,*)RALOCAL(1),RALOCAL(2),RALOCAL(3)                                 

c        print *, indH2O
        if(indH2O) then
          CALL GAGAFE(NATMS, RALOCAL, itagl, FAlocal,
     +                uTotofim, virTot, ETOUT,eQM,natm)
cH2O              this gagafe call is strictly only for water code.
        else
c$$$          CALL GAGAFE(NATMS,RALOCAL,FAlocal,NATMS,RALOCAL,FAlocal,
c$$$     +                uTotofim,virTot,1,1)
        endif

c       print *,uTotofim

c  ----------------

C  INTERACTIONS BETWEEN ATOMS OF TYPE 'ITYPE' WITH THEMSELVES:                            
cH2O      ITPTR=1
cH2O      DO 200 ITYPE=1,NATYPE
cH2O            iatshift1=(ITPTR-1)/3
cH2O            iatshift2=iatshift1
cH2O            CALL GAGAFE(NATMS(ITYPE),RALOCAL(ITPTR),FA(ITPTR),
cH2O     X            NATMS(ITYPE),RALOCAL(ITPTR),FA(ITPTR),
cH2O     X            utotITYPE,virialITYPE,
cH2O     X            ITYPE,1)
cH2O          UT(ITYPE*(ITYPE+1)/2)=utotITYPE
cH2O          virial(ITYPE*(ITYPE+1)/2)=virialITYPE
cH2O        ITPTR=ITPTR+3*NATMS(ITYPE)
cH2O200     CONTINUE
C                                                                                         
c                                                                                         
C  INTERACTIONS BETWEEN ATOMS OF TYPE ITYPE AND JTYPE:                                    
C                                                                                         
cH2O      IF(NATYPE.LE.1) GO TO 500
cH2O      IPOT=NATYPE
cH2O      JTPTR=3*NATMS(1)+1
cH2O      DO 400 JTYPE=2,NATYPE
cH2O        ITPTR=1
cH2O        DO 300 ITYPE=1,JTYPE-1
cH2O          IPOT=IPOT+1
cH2O          iatshift1=(ITPTR-1)/3
cH2O          iatshift2=(JTPTR-1)/3
cH2O          CALL GAGAFE(NATMS(ITYPE),RALOCAL(ITPTR),FA(ITPTR),
cH2O     X            NATMS(JTYPE),RALOCAL(JTPTR),FA(JTPTR),
cH2O     X            utotITYPE,virialITYPE,
cH2O     X            IPOT,0)
cH2O          UT(ITYPE+JTYPE*(JTYPE-1)/2)=utotITYPE
cH2O          virial(ITYPE+JTYPE*(JTYPE-1)/2)=virialITYPE
cH2O          ITPTR=ITPTR+3*NATMS(ITYPE)
cH2O300       CONTINUE
cH2O        JTPTR=JTPTR+3*NATMS(JTYPE)
cH2O400     CONTINUE

cH2O500   CONTINUE

c  ----------------

          fact = 1. / float(nimFPI(1))
          utotim(im) = uTotofim * fact
          do na = 1, natoms
             ind1 = 3 * (na-1)
             ind2 = 3 * ((na-1) * nimFPI(1) + im - 1)
             fa(ind2+1) = falocal(ind1+1) * fact
             fa(ind2+2) = falocal(ind1+2) * fact
             fa(ind2+3) = falocal(ind1+3) * fact
          enddo

cw          write(lunout,*) 'Forces:  utotim(',im,') is ',uTotofim

250     continue
c                close loop over FPI images.

cw    write(lunout,*)'  Forces: after GAGAFE, FA for atom 2 = ',fa(4),fa(5),fa(6)

cPres:                                                                                    
cH2O      virTOT=0.
cH2O      npotint=NATYPE*(NATYPE+1)/2
cH2O      do i=1,npotint
cH2O        virTOT=virTOT+virial(i)
cH2O      enddo

c      UTOT=0.
c      DO 600 I=1,NATYPE*(NATYPE+1)/2
c        UTOT=UTOT+UT(I)
c600     continue

        utot=0.0
        do im = 1, nimFPI(1)
           utot = utot + utotim(im)
        enddo
        if(nFPI .eq. 0) nimFPI(1)=0

c      utot = utotim(im)
c      print *,utotim(im)
c      print *,utotim(1)
c      print *,nimFPI(1)
c      print *, utot

cH2O      totpair=utot

cH2O       CALL EMBED(FRHOTOT,embvir)

cw      write(lunout,*) '  FRHOTOT = ',FRHOTOT,'  embvir = ',embvir                       
cw      write(lunout,*) '  Forces:  after Embed, FA for atom 2 = ',fa(4),fa(5),fa(6)

cPres                                                                                     
cH2O      virTOT = virTOT + embvir
cH2O      UTOT = UTOT + FRHOTOT


c  reset the number of atoms of each type:
      if(nFPI .gt. 0) then
        natoms = natomsold
        do it = 1, natype
          natms(it) = natmsold(it)
        enddo
      endif

c    ------------------------------------------------------------                         
cFPI:  Feynman Path Integrals:  add the potential and force due to interaction            
c                               between images.                                           

       If(nFPI .gt. 0) then
         if (Inudge) then
           call SADDLEFIND
         else
           call FPIforce
         endif
       endif

       if(nFPIcnt .gt. 0) call centroidconstr

c    -------------------------------------------------------------                        


c  Scale force so that Force = actual force / ax :                                        
      TOTFOR = 0.0
      do i = 1, 3 * NATOMS
         TOTFOR = TOTFOR + FA(i)
         FA(i) = FA(i) / ax
      end do

c  Save to output array for python interface
      do i = 1,3*NATOMS
         FAOUT(i) = FA(i)
      end do

c  Write forces and potential for debugging:                                              
cw        write(lunout,*)                                                                 
cw        write(lunout,*)' Forces:    total force = ',TOTFOR                              
cw        write(lunout,*)'   atom:    FAx:     FAy:    FAz:'                              
cw	do ia=1,natoms                                                                         
cw          write(lunout,901) ia,FA(3*(ia-1)+1),FA(3*(ia-1)+2),FA(3*(ia-1)+3)             
cw901       format(' ',i7,3g14.4)                                                         
cw        enddo                                                                           
cw      write(lunout,800) UTOT,(UT(iu),iu=1,NATYPE*(NATYPE+1)/2)                          
cw800   format('     UTOT = ',g20.10/'  UT(ipot) = ',5g12.4)                              
cw      write(lunout,801) totpair,frhotot                                                 
cw801   format('     totpair = ',g12.4,'   frhotot = ',g12.4)                             

c    -------------------------------------------------------------                        
c   Check potential energy per atom by comparing sum with UTOT:                           
cH2O      sum=0.
cH2O      do i=1,natoms
cH2O        Potpat(i)=0.5*Potenpat(i)

cc        write(lunout,*) '   FORCE:   potpat(',i,') = ',Potpat(i)                        

cH2O        sum=sum+Potpat(i)
cH2O      enddo
cp      if(abs(sum-UTOT) .gt. 0.1e-6) write(lunout,802) sum,UTOT                          
cp802   format('  WARNING from FORCE:   sum Potenpat = ',g18.10,                          
cp     +       '  .ne.  UTOT = ',g18.10)                                                  
cc      write(lunout,801) virTOT,(virial(k),k=1,npotint)                                  
cc801   format('   virTOT = ',g20.10,'  VIRIALs:',5g20.10)                                

      RETURN
      END
