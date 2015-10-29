c----------------------------------------------------------------------+
c     Calculate derivatives of the electric field.                     |
c----------------------------------------------------------------------+
      subroutine calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2,
     $     d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)

      implicit none
      include '../commonblks/parameters.cmn'
      integer nM, NC, NCz
      real*8 rCM(3,maxCoo/3), a(3), a2(3)
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)

      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
      real*8 d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3), rMax2 

      real*8 d1d(3), d2d(3,3), d3d(3,3,3)
      real*8 d4d(3,3,3,3), d5d(3,3,3,3,3) 

      real*8 d1a(3), d2a(3,3), d3a(3,3,3)
      real*8 d4a(3,3,3,3), d5a(3,3,3,3,3) 

      real*8 d(3), q(3,3), o(3,3,3), h(3,3,3,3), fsf(3,maxCoo/3)

      integer i, j, k, l, s, n, m, ii, nx, ny, nz
      integer in2(2), in3(3), in4(4), in5(5)
      real*8 re(3), dr(3), r1, r2, swFunc, dSdr
      logical*1 iSlab

ctimming
      integer*8 ti, tf, irtc
      real*8 t1, t2, t3, t4, t5, t6, t7
ctimming

      do n = 1, nM
         do i = 1, 3
            d1v(i,n) = 0.d0
            fsf(i,n) = 0.d0
            do j = i, 3                           
               d2v(i,j,n) = 0.d0
               do k = j, 3                              
                  d3v(i,j,k,n) = 0.d0
                  do l = k, 3
                     d4v(i,j,k,l,n) = 0.d0
                     do s = l, 3
                        d5v(i,j,k,l,s,n) = 0.d0
                     end do
                  end do
               end do
            end do
         end do
      end do

c$$$      t1 = 0.d0  
c$$$      t2 = 0.d0  
c$$$      t3 = 0.d0  
c$$$      t4 = 0.d0  
c$$$      ti = irtc()

      NCz = NC
      if (iSlab) NCz = 0

      do n = 1, nM
         do m = 1, nM
            
            do nx = -NC, NC
               re(1) = a(1) * nx
               do ny = -NC, NC
                  re(2) = a(2) * ny
                  do nz = -NCz, NCz
                     re(3) = a(3) * nz
                  
                     if ( (n.eq.m) .and. (nx.eq.0) .and. (ny.eq.0) .and.
     $                    (nz.eq.0)) goto 11

                     do i = 1, 3
                        dr(i) = rCM(i,n) - rCM(i,m)
                        if (dr(i) .gt. a2(i)) then
                           dr(i) = dr(i) - a(i)
                        else if (dr(i) .lt. -a2(i)) then
                           dr(i) = dr(i) + a(i)
                        end if
                        dr(i) = dr(i) + re(i)
                     end do

                     r2 = dr(1)**2 + dr(2)**2 + dr(3)**2 

                     if (r2 .gt. rMax2) goto 11
                     r1 = sqrt(r2)
                     call SFdsf(r1, swFunc, dSdr)
                     

                     do i = 1, 3
                        d(i) = dpole(i,m)
                     end do
c                     call dDpole(d, dr, d1d, d2d, d3d, d4d, d5d)
                     call dDpole(d, dr, d1a, d2a, d3a, d4a, d5a)


c$$$                     print *, n, m
c$$$                     print *, rCM(2,n), rCM(2,m), a(2)
c$$$
c$$$                     print *, d(1), d(2), d(3)
c$$$                     print *, dr(1), dr(2), dr(3)
c$$$                     print *, d1a(1), d1a(2), d1a(3)
c$$$                     stop

                     do j = 1, 3
                        do i = 1, 3
                           q(i,j) = qpole(i,j,m)
                        end do
                     end do
                     call dQpole(q, dr, d1d, d2d, d3d, d4d, d5d)
                     call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d,
     $                    d3d, d4d, d5d)

c$$$                     print *, d1a(1), d1a(2), d1a(3)
c$$$                     print *, d1d(1), d1d(2), d1d(3)

                     do k = 1, 3
                        do j = 1, 3
                           do i = 1, 3
                              o(i,j,k) = opole(i,j,k,m)
                           end do
                        end do
                     end do
                     call dOpole(o, dr, d1d, d2d, d3d, d4d, d5d)
                     call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d,
     $                    d3d, d4d, d5d)

c$$$                     print *, d1a(1), d1a(2), d1a(3)
c$$$                     print *, d1d(1), d1d(2), d1d(3)

                     do l = 1, 3
                        do k = 1, 3
                           do j = 1, 3
                              do i = 1, 3
                                 h(i,j,k,l) = hpole(i,j,k,l,m)
                              end do
                           end do
                        end do
                     end do
                     call dHpole(h, dr, d1d, d2d, d3d, d4d, d5d)
                     call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d,
     $                    d3d, d4d, d5d)
                     call addDeriv(d1v, d2v, d3v, d4v, d5v, d1a, d2a,
     $                    d3a, d4a, d5a, n, swFunc)

c$$$                     print *, d1a(1), d1a(2), d1a(3)
c$$$                     print *, d1v(1,1), d1v(2,1), d1v(3,1)
c$$$                     stop


                     call addSwitchingForce(d1a, d2a, d3a, d4a, n,
     $                    dSdr, dr, r1, dpole, qpole, opole, hpole, fsf
     $                    )

 11               end do
               end do
            end do
         end do
      end do
c$$$      tf = irtc()
c$$$      t4 = (tf-ti) * 1e-9
c$$$      
c$$$      print '(A,f15.6)', 'Calculation: ', t4

c$$$      print '(A,f15.6)', '      Dipoles: ', t1
c$$$      print '(A,f15.6)', '  Quadrupoles: ', t2
c$$$      print '(A,f15.6)', '    Octopoles: ', t3
c$$$      print '(A,f15.6)', 'hexadecapoles: ', t4


c     Copy all the permutations. (Is this really necessary??)
c$$$      ti = irtc()
      do i = 1, 3
         do j = 1, 3
            in2(1) = i
            in2(2) = j
            call insertIN(in2, 2)
            
            do n = 1, nM
               d2v(i,j,n) = d2v(in2(1), in2(2), n)
            end do

            do k = 1, 3
               do ii = 1, 2
                  in3(ii) = in2(ii)
               end do
               in3(3) = k
               call insertIN(in3, 3)

               do n = 1, nM
                  d3v(i,j,k,n) = d3v(in3(1),in3(2),in3(3),n)
               end do

               do l = 1, 3
                  do ii = 1, 3
                     in4(ii) = in3(ii)
                  end do
                  in4(4) = l
                  call insertIN(in4, 4)

                  do n = 1, nM
                     d4v(i,j,k,l,n) = d4v(in4(1),in4(2),in4(3),in4(4),n)
                  end do
                     
                  do m = 1, 3
                     do ii = 1, 4
                        in5(ii) = in4(ii)
                     end do
                     in5(5) = m
                     call insertIN(in5, 5)
                     
                     do n = 1, nM
                        d5v(i,j,k,l,m,n) = d5v(in5(1),in5(2),in5(3)
     $                       ,in5(4),in5(5),n)
                     end do

                  end do
               end do
            end do
         end do
      end do
c$$$      tf = irtc()
c$$$      t4 = (tf-ti) * 1e-9
c$$$      print '(A,f15.6)', 'Permutations: ', t4
c$$$      stop

      return 
      end
c-----------------------------------------------------------------------
      subroutine insertIN(i, n)

      implicit none
      integer i(*), j, k, iaux, n

      if (n .ge. 2) then
         iaux = i(n)
         j = n-1
         do while (i(j) .gt. iaux .and. j .ge. 1)
            i(j+1) = i(j)
            j = j - 1
         end do
         i(j+1) = iaux
      end if

      return
      end
c-----------------------------------------------------------------------
      function delta(i, j)
      implicit none
      integer i, j, delta

      if (i.eq.j) then
         delta = 1
      else
         delta = 0
      end if

      return
      end
c-----------------------------------------------------------------------
      subroutine addDeriv(d1v, d2v, d3v, d4v, d5v, d1d, d2d, d3d, d4d,
     $     d5d, n, swFunc) 

      implicit none
      include '../commonblks/parameters.cmn'

      real*8 d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
      real*8 d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3) 

      real*8 d1d(3), d2d(3,3), d3d(3,3,3)
      real*8 d4d(3,3,3,3), d5d(3,3,3,3,3) 
      real*8 swFunc
      
      integer n, i, j, k, l, s

      do i = 1, 3
         d1v(i,n) = d1v(i,n) + d1d(i) * swFunc
         do j = i, 3                           
            d2v(i,j,n) = d2v(i,j,n) + d2d(i,j) * swFunc
            do k = j, 3                              
               d3v(i,j,k,n) = d3v(i,j,k,n) + d3d(i,j,k) * swFunc
               do l = k, 3
                  d4v(i,j,k,l,n) = d4v(i,j,k,l,n) + d4d(i,j,k,l)
     $                 * swFunc
                  do s = l, 3
                     d5v(i,j,k,l,s,n) = d5v(i,j,k,l,s,n) + d5d(i,j,k,l,s
     $                    ) * swFunc
                  end do
               end do
            end do
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      subroutine addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d, d3d, d4d,
     $     d5d) 

      implicit none
      include '../commonblks/parameters.cmn'

      real*8 d1a(3), d2a(3,3), d3a(3,3,3)
      real*8 d4a(3,3,3,3), d5a(3,3,3,3,3) 

      real*8 d1d(3), d2d(3,3), d3d(3,3,3)
      real*8 d4d(3,3,3,3), d5d(3,3,3,3,3) 

      integer n, i, j, k, l, s

      do i = 1, 3
         d1a(i) = d1a(i) + d1d(i)
         do j = i, 3                           
            d2a(i,j) = d2a(i,j) + d2d(i,j)
            do k = j, 3                              
               d3a(i,j,k) = d3a(i,j,k) + d3d(i,j,k)
               do l = k, 3
                  d4a(i,j,k,l) = d4a(i,j,k,l) + d4d(i,j,k,l)
                  do s = l, 3
                     d5a(i,j,k,l,s) = d5a(i,j,k,l,s) + d5d(i,j,k,l,s)
                  end do
               end do
            end do
         end do
      end do

      return
      end
c-----------------------------------------------------------------------
      subroutine addSwitchingForce(d1a, d2a, d3a, d4a, n, dSdr, dr,
     $     r1, dpole, qpole, opole, hpole, fsf)

      implicit none
      include '../commonblks/parameters.cmn'

      real*8 d1a(3), d2a(3,3), d3a(3,3,3), d4a(3,3,3,3)
      real*8 dSdr, dr(3), r1, fsf(3,maxCoo/3), u
      real*8 dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
      real*8 opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
      
      integer i, ii, j, k, l, in2(2), in3(3), in4(4), n

c     Copy all the permutations. 
      do i = 1, 3
         do j = 1, 3
            in2(1) = i
            in2(2) = j
            call insertIN(in2, 2)
            
            d2a(i,j) = d2a(in2(1), in2(2))

            do k = 1, 3
               do ii = 1, 2
                  in3(ii) = in2(ii)
               end do
               in3(3) = k
               call insertIN(in3, 3)

               d3a(i,j,k) = d3a(in3(1),in3(2),in3(3))

               do l = 1, 3
                  do ii = 1, 3
                     in4(ii) = in3(ii)
                  end do
                  in4(4) = l
                  call insertIN(in4, 4)

                  d4a(i,j,k,l) = d4a(in4(1),in4(2),in4(3),in4(4))
                     
               end do
            end do
         end do
      end do
      
      u = 0.d0
      do i = 1, 3
         u = u + d1a(i) * dpole(i,n)
         do j = 1, 3
            u = u + d2a(i,j) * qpole(i,j,n) / 3.d0
            do k = 1, 3
               u = u + d3a(i,j,k) * opole(i,j,k,n) / 15.d0
               do l = 1, 3
                  u = u + d4a(i,j,k,l) * hpole(i,j,k,l,n) / 105.d0
               end do
            end do
         end do
      end do

      u = -u * dSdr / r1
      do i = 1, 3
         fsf(i,n) = fsf(i,n) + u * dr(i)
      end do

      return
      end
