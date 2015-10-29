C AA 29-06-2004
C Read in parameters for the potential

      SUBROUTINE potinit()
        IMPLICIT NONE
        INTEGER i,ip,stat
        CHARACTER(LEN=80) text

        INCLUDE '../commonblks/potparam.cmn'
        INCLUDE '../commonblks/compotent.cmn'           

        OPEN(421,FILE='param.pot',STATUS='old',ACTION='read')
        READ(421,*) (fsti(i),text , i=1,nl)        

        READ(421,'(/)')
        READ(421,*) (rcut(i),i=1,3)
        READ(421,*) (rskin(i),i=1,3)
        READ(421,*) indf1

        npotpar=0
        DO ip=1,MAXPOTPAR
           READ(421,*,IOSTAT=stat) (potpar(ip,i),i=1,3)
           IF (stat /= 0) EXIT
           npotpar=npotpar+1
        END DO
c        WRITE(*,*) ' FROM POTINIT ... npotpar = ',npotpar

        CLOSE(421)

        DO  i=1,3
          rcut2(i)=rcut(i)**2
          rskin2(i)=rskin(i)**2
        END DO

      END 
