c-----------------------------------------------------------------------
      subroutine molForce3(f1, f2, f3, r1, r2, r3, fTot, tt, flag)

      implicit none
      integer i, j, ii, jj, ll
      real*8 r1(3), r2(3), r3(3), rc1(3), rc2(3), rc3(3) 
      real*8 f1(3), f2(3), f3(3), tor(3), fidd(3)

      real*8 b(6), c(6), y(6,6)
      real*8 tt(3), ftot(3), t1(3), t2(3), t3(3)
      integer flag
      real*8 f1a(3), f2a(3), f3a(3)
c      data mF / 0.352d0, -0.864d0, -0.36d0, 0.36d0, 0.48d0, -0.8d0,
c     $     .864d0, 0.152d0, 0.48d0 /
c      data mB / 0.352d0, 0.36d0, 0.864d0, -0.864d0, 0.48d0, 0.152d0,
c     $     -0.36d0, -0.8d0, 0.48d0/

      real*8 mF(3,3,3), mB(3,3,3), sq2, oh
      parameter (sq2=0.707106781186547524d0, oh=0.5d0)
c$$$      data mF / -sq2, oh, -oh, 0.d0, sq2, sq2, sq2, oh, -oh,
c$$$     $     -oh, -sq2, -oh, oh, -sq2, oh, -sq2, 0.d0, sq2,
c$$$     $     sq2, sq2, 0.d0, oh, -oh, -sq2, -oh, oh, -sq2 /
c$$$
c$$$      data mB / -sq2, 0.d0, sq2, oh, sq2, oh, -oh, sq2, -oh,
c$$$     $     -oh, oh, -sq2, -sq2, -sq2, 0.d0, -oh, oh, sq2,
c$$$     $     sq2, oh, -oh, sq2, -oh, oh, 0.d0, -sq2, -sq2 /

      data mF / -0.707106781186547524d0, 0.5d0, -0.5d0, 0.d0, 0
     $     .707106781186547524d0, 0.707106781186547524d0, 0
     $     .707106781186547524d0, 0.5d0, -0.5d0,-0.5d0, -0
     $     .707106781186547524d0, -0.5d0, 0.5d0, -0.707106781186547524d0
     $     , 0.5d0, -0.707106781186547524d0, 0.d0, 0
     $     .707106781186547524d0,0.707106781186547524d0, 0
     $     .707106781186547524d0, 0.d0, 0.5d0, -0.5d0, -0
     $     .707106781186547524d0, -0.5d0, 0.5d0, -0.707106781186547524d0
     $     /

      data mB / -0.707106781186547524d0, 0.d0, 0.707106781186547524d0, 0
     $     .5d0, 0.707106781186547524d0, 0.5d0, -0.5d0, 0
     $     .707106781186547524d0, -0.5d0,-0.5d0, 0.5d0, -0
     $     .707106781186547524d0, -0.707106781186547524d0, -0
     $     .707106781186547524d0, 0.d0, -0.5d0, 0.5d0, 0
     $     .707106781186547524d0,0.707106781186547524d0, 0.5d0, -0.5d0,
     $     0.707106781186547524d0, -0.5d0, 0.5d0, 0.d0, -0
     $     .707106781186547524d0, -0.707106781186547524d0/

      call inv6(r1, r2, r3, y, flag)

c$$$c Total force

c$$$      if (flag .gt. 0) then
c$$$      do i = 1, 3
c$$$         fTot(i) = f1(i) + f2(i) + f3(i) + fidd(i)
c$$$      end do
c$$$      call cross (rc1, f1, t1)
c$$$      call cross (rc2, f2, t2)
c$$$      call cross (rc3, f3, t3)
c$$$c Total torque
c$$$      do i = 1, 3
c$$$         tt(i) = t1(i) + t2(i) + t3(i) + tor(i)
c$$$      end do
c$$$
c$$$c----- Start Debugging ----
c$$$c      print *, 'NEW --------'
c$$$      print *
c$$$      print '(6f19.14)', tt(1), tt(2), tt(3), ftot(1), ftot(2), 
c$$$     $     ftot(3)
c$$$      print '(9g15.5)', (f1(ll), ll=1,3), (f2(ll), ll=1,3), 
c$$$     $     (f3(ll), ll=1,3) 
c$$$      end if
c----- End Debugging ----

c     The determinant was zero, therefore, we rotated the space. We have
c     to rotate also the force and torque and then bring things back to
c     the original orientation.
      if (flag .gt. 0) then
         do ii = 1, 3
            f1a(ii) = 0.d0
            f2a(ii) = 0.d0
            do jj = 1, 3
               f1a(ii) = f1a(ii) + mF(ii,jj,flag) * fTot(jj)
               f2a(ii) = f2a(ii) + mF(ii,jj,flag) * tt(jj)
            end do
         end do
         do ii = 1, 3
            ftot(ii) = f1a(ii)
            tt(ii) = f2a(ii)
         end do
      end if

      do i = 1, 3
         b(i) = fTot(i)
         b(i+3) = tt(i)
      end do
            
      do i = 1, 6
         c(i) = 0.d0
         do j = 1, 6
            c(i) = c(i) + y(i,j) * b(j)
         end do
      end do

c     New forces...

      f1(1) = c(1)
      f1(2) = 0.d0
      f1(3) = c(3)

      f2(1) = 0.d0
      f2(2) = c(2)
      f2(3) = c(6)

      f3(1) = c(4)
      f3(2) = c(5)
      f3(3) = 0.d0

c     Here we bring the forces to the original orientation of the space. 
      if (flag .gt. 0) then
         do ii = 1, 3
            f1a(ii) = 0.d0
            f2a(ii) = 0.d0
            f3a(ii) = 0.d0
            do jj = 1, 3
               f1a(ii) = f1a(ii) + mB(ii,jj,flag) * f1(jj)
               f2a(ii) = f2a(ii) + mB(ii,jj,flag) * f2(jj)
               f3a(ii) = f3a(ii) + mB(ii,jj,flag) * f3(jj)
            end do
         end do
         do ii = 1, 3
            f1(ii) = f1a(ii)
            f2(ii) = f2a(ii)
            f3(ii) = f3a(ii)
         end do
      end if

c----- Start Debugging ----
c$$$      if (flag .gt. 0) then
c$$$      call cross (r1, f1, t1)
c$$$      call cross (r2, f2, t2)
c$$$      call cross (r3, f3, t3)
c$$$c Total torque
c$$$      do i = 1, 3
c$$$         tt(i) = t1(i) + t2(i) + t3(i)
c$$$      end do
c$$$c Total force
c$$$      do i = 1, 3
c$$$         ftot(i) = f1(i) + f2(i) + f3(i)
c$$$      end do
c$$$      print '(6f19.10)', tt(1), tt(2), tt(3), ftot(1), ftot(2), 
c$$$     $     ftot(3)
c$$$c      print '(9g15.5)', (f1(ll), ll=1,3), (f2(ll), ll=1,3), 
c$$$c     $     (f3(ll), ll=1,3) 
c$$$      print *, '------------'
c$$$      end if
c----- End Debugging ----

      return
      end
