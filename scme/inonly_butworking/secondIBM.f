c  This routine needs to be taylored for the machine this code is being run on:

       SUBROUTINE SECOND(FLO)
c           (on the CRAY the system timer is called SECOND so this definition
c            of the function should then be commented out.) 

       implicit real*8 (a-h,o-z)
       real*4 secvax,floc,timeunix(2),tzero

c           tzero=0.0
c           floc=cputime(tzero)
c           flo=dble(floc)
c                for the CONVEX using the VECLIB library.

cUNIX          flo=etime(timeunix)
c                for UNIX.

c        call clockx(flo)
c        flo=flo*1.0d-06
c                for the IBM 3090.

cvax       secvax=0.0
cvax       floc=secnds(secvax)
cvax       flo=dble(floc)
c                for the microVAX 3600, it gives REAL time not CPU time.
                  
       return
       end
 
 
