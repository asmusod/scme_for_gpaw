c  Doesn't include piston motion.  So, only constant volume.

c       QUICKMIN, modifies the velocity of the atoms so as to
c       minimize the energy of a config quickly (just a local
c       minimization) by zeroing the component of the velocity perp. to
c       the force vector, and also zeroing the parallel
c       component if it's in the opposite direction to the force vector.

      subroutine QUICKMIN

      implicit real*8 (a-h,o-z)

      include '../commonblks/parameters.cmn'
      include '../commonblks/comgeom.cmn'
      include '../commonblks/comconf.cmn'
      include '../commonblks/comtime.cmn'
      include '../commonblks/combaths.cmn'
      include '../commonblks/comluns.cmn'
      include '../commonblks/comoutcntrl.cmn'

c     Calculate the dot product of the 3N-dimensional velocity and force
c     vectors. (N is the # of atoms, the quickmin routine does all of
c     them.) 

      vdotf = 0.
      do icoo = 1,3*NATOMS
         vdotf = vdotf + va(icoo)*fa(icoo)
      enddo

      if (vdotf .le. 0.) then
c     If velocity is in wrong direction, zero it.
c$$$         write(lunqckmin,241)
 241     format('   zero velocities, vdotf .le. 0.')
         fmag2 = 0.
         fmax = 0.
         fmin = 0.
         do icoo = 1,3*NATOMS
            va(icoo) = 0.

c     find maximum force, just to monitor the convergence:
            fmag2 = fmag2 + fa(icoo)**2
            if(fa(icoo) .lt. fmin) then
               fmin=fa(icoo)
               icoomin=icoo
            endif
            if(fa(icoo) .gt. fmax) then
               fmax=fa(icoo)
               icoomax=icoo
            endif
         enddo
         
      else
c     Project out the parallel component of velocity; this becomes new
c     velocity. 

         fmag2 = 0.
         fmax = 0.
         fmin = 0.
         do icoo = 1,3*NATOMS
            fmag2 = fmag2 + fa(icoo)**2
            if(fa(icoo) .lt. fmin) then
               fmin=fa(icoo)
               icoomin=icoo
            endif
            if(fa(icoo) .gt. fmax) then
               fmax=fa(icoo)
               icoomax=icoo
            endif
         enddo
cc        fmag = sqrt(fmag2)
c     This is length of force vector, 
c     vdotf/fmag2 is length of new velocity vector.

         if (fmag2 .gt. 0.) then
            scale = vdotf/fmag2
            do icoo = 1,3*NATOMS
               va(icoo) = fa(icoo)*scale
            enddo
         endif

      endif
c                close `if (vdotf .le. 0.) then'


cw         write(lunqckmin,*)' QUICKMIN: time=',time,' fmag2=',fmag2
cw         write(lunqckmin,*)'   fa(4..5..6) = ',fa(4),fa(5),fa(6)
cw           write(lunout,233) 1,(va(3*k),k=1,natms(1))
cw           write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))
cw233        format('   new velocities for type',i1,'atoms:'/
cw     +                5g12.6/5g12.6/5g12.6/5g12.6)

      iatmax=icoomax/3+1
      if(-fmin .gt. fmax) then
         iatmax=icoomin/3+1
         fmax=-fmin
      endif
      ind=3*iatmax
c$$$      IF(Jstatwr+1 .ge. Nstatwr) then
c     only write the info into the outqckmin.dat file if its
c     time to write info on trajectory
c$$$         write(lunqckmin,240) 
c$$$     +        time,sqrt(fmag2/Natoms)*ax,fmax*ax,iatmax+1,
c$$$     +        Ra(ind+1)*ax,Ra(ind+2)*ax,Ra(ind+3)*ax,
c$$$     +        Va(ind+1)*ax,Va(ind+2)*ax,Va(ind+3)*ax,
c$$$     +        Fa(ind+1)*ax,Fa(ind+2)*ax,Fa(ind+3)*ax
c$$$      endif
 240  format(' ',f10.5,' ',2g12.4,i5,'   ',3f9.4,' ',3f9.4,' ',
     +     3f11.7)

c     Check to see if convergence is reached:
      if(fmax .lt. tolforce) Iqkmconv=.true.

      return
      end
