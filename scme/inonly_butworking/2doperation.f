c   md2d operation:                                                                       
c      this subroutine zeros the velocity and force in z direction for 2D simul.          

       subroutine md2dop(v2md2d)

       implicit real*8 (a-h,o-z)

       include '../commonblks/parameters.cmn'
       include '../commonblks/comgeom.cmn'
       include '../commonblks/comconf.cmn'
       include '../commonblks/comtime.cmn'
       include '../commonblks/comluns.cmn'

cd      write(6,*) '  MD2D:  time = ',time,'   natoms = ',natoms                          

c     Zero both the velocity and force on all atoms in the z direction:                   
         do 136 i=1,natoms
            v2md2d=v2md2d+va(3*(i-1)+3)**2
c                'v2md2d' keeps track of the total kinetic energy taken                   
c                         from the atoms.  Just for curiosity.                            
            va(3*(i-1)+3)=0.0
            fa(3*(i-1)+3)=0.0
136         continue

cd         write(lunout,*)  '  MD2D:  natoms = ',natoms                                   
cd           write(lunout,*) '  After md2d operation: time=',time                         
cd           write(lunout,233) 1,(va(3*k),k=1,natms(1))                                   
cd           write(lunout,233) 2,(va(3*(natms(1)+k)),k=1,natms(2))                        
cd233        format('   z component of velocities for type',i1,'atoms:'/                  
cd     +                5g15.6/5g15.6/5g15.6/5g15.6)                                      

        return
      end
