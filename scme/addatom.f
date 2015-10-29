c  Addatom:
c      this subroutine adds a new atom to all the lists of atoms.
c      It will become the first atom of its type in the lists (i.e. atom #1
c      if it is of type A, atom # natms(1)+1 if it is of type B)
c
      subroutine addatom(itype,xcoo,ycoo,zcoo,xvel,yvel,zvel,navecoo)
c
       implicit real*8 (a-h,o-z)

       include '../commonblks/parameters.cmn'
       include '../commonblks/comconf.cmn'
       include '../commonblks/comgeom.cmn'
       include '../commonblks/comtime.cmn'
       include '../commonblks/comluns.cmn'
       include '../commonblks/comintlis.cmn'



c    Find where in the lists this atom lends:  (assume this type already exists)
       ibeg=1
c             this is correct if the new atom is of type A.
       do it=1,NATYPE
         if(itype .eq. it) go to 50
         ibeg=ibeg+natms(it)
       enddo
50     continue
c             ibeg is now the atom number of the new atom.


c    Shift the lists of coordinates, velocities and forces to make room for
c    the new atom, start from the high end:
c
             do 195 iresp=ibeg,natoms
                ires=natoms-(iresp-ibeg)

                ra(3*(ires+1)-2)=ra(3*ires-2)
                ra(3*(ires+1)-1)=ra(3*ires-1)
                ra(3*(ires+1))  =ra(3*ires)
c
                va(3*(ires+1)-2)=va(3*ires-2)
                va(3*(ires+1)-1)=va(3*ires-1)
                va(3*(ires+1))  =va(3*ires)
c
                fa(3*(ires+1)-2)=fa(3*ires-2)
                fa(3*(ires+1)-1)=fa(3*ires-1)
                fa(3*(ires+1))  =fa(3*ires)
c
                avera(3*(ires+1)-2)=avera(3*ires-2)
                avera(3*(ires+1)-1)=avera(3*ires-1)
                avera(3*(ires+1))  =avera(3*ires)
c
                sumra(3*(ires+1)-2)=sumra(3*ires-2)
                sumra(3*(ires+1)-1)=sumra(3*ires-1)
                sumra(3*(ires+1))  =sumra(3*ires)
c
                raNoPBC(3*(ires+1)-2)=raNoPBC(3*ires-2)
                raNoPBC(3*(ires+1)-1)=raNoPBC(3*ires-1)
                raNoPBC(3*(ires+1))  =raNoPBC(3*ires)

                itag(ires+1)=itag(ires)
                Potpat(ires+1)=Potpat(ires)


195             continue
c
           natoms=natoms+1
           natms(itype)=natms(itype)+1
c
c    Now add the new atom to the lists:
                ra(3*ibeg-2)=xcoo
                ra(3*ibeg-1)=ycoo
                ra(3*ibeg)  =zcoo
c
                va(3*ibeg-2)=xvel
                va(3*ibeg-1)=yvel
                va(3*ibeg)  =zvel
c
                fa(3*ibeg-2)=0.0
                fa(3*ibeg-1)=0.0
                fa(3*ibeg)  =0.0
c
                avera(3*ibeg-2)=xcoo
                avera(3*ibeg-1)=ycoo
                avera(3*ibeg)  =zcoo
c
                sumra(3*ibeg-2)=xcoo*navecoo
                sumra(3*ibeg-1)=ycoo*navecoo
                sumra(3*ibeg)  =zcoo*navecoo

                raNoPBC(3*ibeg-2)=xcoo
                raNoPBC(3*ibeg-1)=ycoo
                raNoPBC(3*ibeg)  =zcoo

c             this needs to be multiplied by the number of terms in sum so far.
c
                itag(ibeg)=0
c                      tag this atom as an ordinary fully dynamical atom.
c
c       Fix the list of permafrost atoms, now that indices have changed:
          do i=1,naperm
	     if(iperm(i) .ge. ibeg) iperm(i)=iperm(i)+1
	  enddo
c               it is assumed that the new atom does not land in the permafrost.
c       Fix the list of FPI atoms, now that indices have changed:
          do i=1,nFPI
	     if(iFPI(i) .ge. ibeg) iFPI(i)=iFPI(i)+1
	  enddo
c       Fix pointers to atoms:
          if(indzbott .ge. ibeg) indzbott=indzbott+1
          if(indzsurf .ge. ibeg) indzsurf=indzsurf+1
          if(indatwr .ge. ibeg) indatwr=indatwr+1
c
c         make sure neighborlist gets updated:
          indupd=1	  
c                 ------------------------------
c
cd             write(lunout,590) itype,natoms,natms(1),time,
cd     +                         xcoo,ycoo,zcoo,zvel
cd590          format(' Addatom subr:  type',i1,' Natoms =',i5,' Na =',i5,
cd     +              ' time =',f9.3/
cd     +              '   x =',g12.4,' y =',g12.4,' z =',g12.4,
cd     +              '   V0z =',g12.4)
c
         return
	 end
