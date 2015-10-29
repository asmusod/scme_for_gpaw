c   version f:  Add possibility that all atoms are in chains,  indallinchains=1
c   version e:   Hong adds PBC on the FPI force
c   version c:   make spring constant only type dependent.
c                introduce circular FPI chains, with or without centroid
c                constraint. 
c   EAM version b91: add FPI.   calculate potential and force between adjacent
c                               images in FPI and add to FA, Potperat, UTOT,

c     Potential energy between adjacent images of FPI chains is
c                        U = 0.5 * springconst * delr**2 
c     The force is       F = - springconst * delr
c-----------------------------------------------------------------------
      subroutine FPIforce

      implicit real*8 (a-h,o-z)

      dimension del(3),PL(3)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comenergy.cmn'
      
      PL(1)=ax
      PL(2)=ay
      PL(3)=az
c-----------------------------------------------------------------------
c     Loop over FPI chains:
                                                                 
      do 1 ipf=1,nFPI

c     get the type of the atoms in this chain:
         nsofar=0
         do itype=1,natype
            nsofar = nsofar + NATMS(itype)
            if(iFPI(ipf) .lt. nsofar) go to 400
         enddo
 400     continue

         do 2 im = iFPI(ipf), iFPI(ipf) + nimFPI(ipf) - 2
c     loop over images in FPI chain `ipf' skipping the last image:
cPBC:
            do ii=1,3
               del(ii) = (RA(3*(im-1)+ii) - RA(3*(im)+ii)) * ax
               if(abs(del(ii)) .gt. 0.5*PL(ii)) 
     &              del(ii) = del(ii) * (1. - PL(ii)/abs(del(ii))) 
            enddo
            r2 = del(1)*del(1) + del(2)*del(2) + del(3)*del(3)
            utot = utot + 0.5d0 * sprconFPI(itype) * r2

c     add to force on first image:
            FA(3*(im-1)+1) = FA(3*(im-1)+1) - sprconFPI(itype) * del(1)
            FA(3*(im-1)+2) = FA(3*(im-1)+2) - sprconFPI(itype) * del(2)
            FA(3*(im-1)+3) = FA(3*(im-1)+3) - sprconFPI(itype) * del(3)
c     add to force on second image:
            FA(3*(im)+1) = FA(3*(im)+1) - sprconFPI(itype) * (-del(1))
            FA(3*(im)+2) = FA(3*(im)+2) - sprconFPI(itype) * (-del(2))
            FA(3*(im)+3) = FA(3*(im)+3) - sprconFPI(itype) * (-del(3))
 2       continue

c     CircularFPI:
         if((itagFPI(ipf) .ge. 200 .and. itagFPI(ipf) .lt. 400) .or.
     +        (indallinchains .eq. 1 .and. indallfixend .eq. 0)) then
c AllinCh. 
c Here the FPI is circular, add spring interaction between first and
c last image in the FPI:

            im1 = iFPI(ipf)
            im2 = iFPI(ipf) + nimFPI(ipf) - 1

cPBC:
            do ii = 1, 3
               del(ii) = (RA(3*(im1-1) + ii) - RA(3*(im2-1)+ii)) * ax
               if(abs(del(ii)) .gt. 0.5*PL(ii))
     &              del(ii) = del(ii) * (1. - PL(ii)/abs(del(ii)))
            enddo
            r2=del(1)*del(1) + del(2)*del(2) + del(3)*del(3)
            utot = utot + 0.5 * sprconFPI(itype) * r2

c         add to force on first image:
            FA(3*(im1-1)+1) = FA(3*(im1-1)+1) - sprconFPI(itype)*del(1)
            FA(3*(im1-1)+2) = FA(3*(im1-1)+2) - sprconFPI(itype)*del(2)
            FA(3*(im1-1)+3) = FA(3*(im1-1)+3) - sprconFPI(itype)*del(3)
            
c     add to force on second image:
            FA(3*(im2-1)+1) = FA(3*(im2-1)+1)-sprconFPI(itype)*(-del(1))
            FA(3*(im2-1)+2) = FA(3*(im2-1)+2)-sprconFPI(itype)*(-del(2))
            FA(3*(im2-1)+3) = FA(3*(im2-1)+3)-sprconFPI(itype)*(-del(3))
         endif
         
 1    continue
      
      return
      end
 
