c   removeatom:   
c        This subroutine removes the atom number indat from all the vectors
c        that list atoms and does the shifts in the vectors.
c
c
      subroutine removeatom(indat)
c
      implicit real*8 (a-h,o-z)


       include '../commonblks/parameters.cmn'
       include '../commonblks/comconf.cmn'
       include '../commonblks/comgeom.cmn'
       include '../commonblks/comtime.cmn'
       include '../commonblks/comluns.cmn'
       include '../commonblks/comintlis.cmn'


      dimension rdispl(MAXCOO)
 
       if(indat .gt. Natoms) then
          write(6,*) ' ERROR in removeatom:  indat > Natoms, indat = ',indat
	  stop
       endif
       if(indat .eq. Natoms) go to 10

             do 115 ires=indat+1,natoms
c                from comconf.cmn
                RA(3*(ires-1)-2)=RA(3*ires-2)
                RA(3*(ires-1)-1)=RA(3*ires-1)
                RA(3*(ires-1))  =RA(3*ires)
c
                VA(3*(ires-1)-2)=VA(3*ires-2)
                VA(3*(ires-1)-1)=VA(3*ires-1)
                VA(3*(ires-1))  =VA(3*ires)
c
                fa(3*(ires-1)-2)=fa(3*ires-2)
                fa(3*(ires-1)-1)=fa(3*ires-1)
                fa(3*(ires-1))  =fa(3*ires)
c
                avera(3*(ires-1)-2)=avera(3*ires-2)
                avera(3*(ires-1)-1)=avera(3*ires-1)
                avera(3*(ires-1))  =avera(3*ires)
c
                sumra(3*(ires-1)-2)=sumra(3*ires-2)
                sumra(3*(ires-1)-1)=sumra(3*ires-1)
                sumra(3*(ires-1))  =sumra(3*ires)
c
                raNoPBC(3*(ires-1)-2)=raNoPBC(3*ires-2)
                raNoPBC(3*(ires-1)-1)=raNoPBC(3*ires-1)
                raNoPBC(3*(ires-1))  =raNoPBC(3*ires)

c         lists of atoms:
                itag(ires-1)=itag(ires)
                Potpat(ires-1)=Potpat(ires)

115             continue
c
10           continue

              isitperm=0
              do 117 ip=1,naperm
		 if(iperm(ip) .ne. indat) go to 117
c                           so atom 'indat' is a permafrost atom.
                    isitperm=1
                    ipindat=ip      
117              continue
              if(isitperm .eq. 1) then
	        do iper=ipindat+1,naperm
		   iperm(iper-1)=iperm(iper)
                enddo		
                naperm=naperm-1
              endif
cc              naperm=naperm-1
c         shift the list of permafrost atoms:
              do ip=1,naperm
                 if(iperm(ip) .gt. indat) iperm(ip)=iperm(ip)-1
              enddo
c         shift the list of permafrost atoms:
              do ip=1,nFPI
                 if(iFPI(ip) .gt. indat) iFPI(ip)=iFPI(ip)-1
              enddo
c         Fix pointers to atoms:
          if(indzbott .gt. indat) indzbott=indzbott-1
          if(indzsurf .gt. indat) indzsurf=indzsurf-1

c          Neighborlist must be updated the next time GAGAFE is called:
              indupd=1
c
       natoms=natoms-1
c   Find which type the atom is:
       nsum=natms(1)
       do itype=1,NATYPE
         if(indat .le. nsum) go to 50
         nsum=nsum+natms(itype+1)
       enddo
50     continue

       natms(itype)=natms(itype)-1
       if(natms(itype).eq.0) then
c         set some variables to flag that an entire component has been
c          removed
         istpcomp = 1
         istpno = itype
         do 500 ia = itype+1,NATYPE
           natms(ia-1)=natms(ia)
           AMASS(ia-1)=AMASS(ia)
           do 490 ib = 1,80
             CPHEAD(ia-1,ib)=CPHEAD(ia,ib)
490        continue
500      continue
         NATYPE=NATYPE-1
       endif
c
cd       write(6,592) indat,itype,natoms,natms(1),natms(2)
cd592    format('  Removeatoms:  remove atom ',i5,' of type',i2,
cd     +         ' natoms = ',i5,' natms(1) = ',i5,' natms(2)=',i5)
c
       return 
       end

