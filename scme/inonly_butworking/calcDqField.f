c----------------------------------------------------------------------+
c     Derivatives of the quadrupole field                              |
c----------------------------------------------------------------------+
      subroutine dQpole(q, r, d1q, d2q, d3q, d4q, d5q)

      implicit none
      real*8 q(3,3), r(3), d1q(3), d2q(3,3), d3q(3,3,3)
      real*8 d4q(3,3,3,3), d5q(3,3,3,3,3) 

      integer delta, dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
      integer i, j, k, l, s
      real*8 r2, r3, r5, r7, r9, r11, r13, r15, rrq
      real*8 t1, t2, t3, y1, y2, y3, z1, z2, z3, z4, w1, w2, w3, w4
      real*8 v(3)

      real*8 dt2k, dt3k, ddt2kl, dy1l, dt2l, ddt3kl, dy2l, dt3l, dy3l
      real*8 dddt3kls, dt2s, ddt2ls, ddy2ls, ddt2ks, dy1s, dz2s
      real*8 ddt3ks, dy2s, ddt3ls, dt3s, ddy3ls, dz3s, dy3s, dz4s

      r2 = r(1)**2 + r(2)**2 + r(3)**2
      r3 = sqrt(r2) * r2
      r5 = r3 * r2
      r7 = r5 * r2
      r9 = r7 * r2
      r11 = r9 * r2
      r13 = r11 * r2
      r15 = r13 * r2

      r5  =        1.d0 / r5
      r7  =       -5.d0 / r7
      r9  =       35.d0 / r9
      r11 =     -315.d0 / r11
      r13 =     3465.d0 / r13
      r15 =   -45045.d0 / r15

      rrq = 0.d0
      do j = 1, 3
         v(j) = 0.d0
         do i = 1, 3
            rrq    = rrq + q(i,j) * r(i) * r(j)
            v(j)   = v(j) + q(i,j) * r(i)
         end do
      end do

      do i = 1, 3
         d1q(i) = (2.d0 * v(i)) * r5 + (rrq * r(i)) * r7

         do j = i, 3
            dij = delta(i,j)

c     2nd derivative
            t1 = 2.d0 * q(i,j)
            t2 = 2.d0 * (v(i)*r(j) + v(j)*r(i)) + rrq * dij
            t3 = rrq * r(i)*r(j)

            d2q(i,j) = t1 * r5 + t2 * r7 + t3 * r9

            do k = j, 3
               dik = delta(i,k)
               djk = delta(j,k)
               
c     3rd derivative
               dt2k = 2.d0 * (v(i)*djk + v(j)*dik + v(k) * dij) + 2.d0
     $              * (q(i,k)*r(j) + q(j,k)*r(i))
               dt3k = 2.d0 * v(k)*r(i)*r(j) + rrq * (r(i)*djk + r(j)*dik
     $              )
               y1 = t1*r(k) + dt2k
               y2 = t2*r(k) + dt3k
               y3 = t3*r(k)
               
               d3q(i,j,k) = y1*r7 + y2*r9 + y3*r11

               do l = k, 3
                  dil = delta(i,l)
                  djl = delta(j,l)
                  dkl = delta(k,l)
               
c     4th derivative
                  ddt2kl = 2.d0 * (q(i,l)*djk + q(j,l)*dik + q(k,l)
     $                 * dij) + 2.d0 * (q(i,k)*djl + q(j,k)*dil)
                  dy1l = t1*dkl + ddt2kl

                  z1 = dy1l

                  dt2l = 2.d0 * (v(i)*djl + v(j)*dil + v(l) * dij) + 2
     $                 .d0* (q(i,l)*r(j) + q(j,l)*r(i))


                  ddt3kl = 2.d0 * q(k,l)*r(i)*r(j) + 2.d0 * v(k)* (r(i)
     $                 *djl + dil*r(j) )
     $                 + 2.d0 * v(l) * (r(i)*djk + r(j)*dik)
     $                 + rrq * (dil*djk + djl*dik)
               
                  dy2l = dt2l*r(k) + t2*dkl + ddt3kl

                  z2 = y1*r(l) + dy2l
                  

                  dt3l = 2.d0 * v(l)*r(i)*r(j) + rrq * (r(i)*djl + r(j)
     $                 *dil)

                  dy3l = dt3l*r(k) + t3 * dkl
                  z3 = y2*r(l) + dy3l

                  z4 = y3*r(l)
                  d4q(i,j,k,l) =  z1*r7 + z2*r9 + z3*r11 + z4*r13
                  
                  do s = l, 3
                     dis = delta(i,s)
                     djs = delta(j,s)
                     dks = delta(k,s)
                     dls = delta(l,s)


c     5th derivative
                     
                     dddt3kls = 2.d0 * q(k,l)*(r(i)*djs + dis*r(j))
     $                    + 2.d0 * q(k,s) * (r(i)*djl + dil*r(j))
     $                    + 2.d0 * v(k) * (dis*djl + dil*djs)
     $                    + 2.d0 * q(l,s) * (r(i)*djk + r(j)*dik)
     $                    + 2.d0 * v(l) * (dis*djk + djs*dik)
     $                    + 2.d0 * v(s) * (dil*djk + djl*dik)

                     dt2s = 2.d0 * (v(i)*djs + v(j)*dis + v(s) * dij) +
     $                    2.d0* (q(i,s)*r(j) + q(j,s)*r(i))

                     ddt2ls = 2.d0 * (q(i,s)*djl + q(j,s)*dil + q(l,s)
     $                    * dij) + 2.d0* (q(i,l)*djs + q(j,l)*dis)

                     ddy2ls = ddt2ls*r(k) + dt2l*dks + dt2s*dkl +
     $                    dddt3kls  

                     ddt2ks = 2.d0 * (q(i,s)*djk + q(j,s)*dik + q(k,s)
     $                    * dij) + 2.d0* (q(i,k)*djs + q(j,k)*dis)

                     dy1s = t1*dks + ddt2ks

                     dz2s = y1*dls + dy1s*r(l) + ddy2ls
                     w1 = z1*r(s) + dz2s
c----------
                     ddt3ks = 2.d0 * q(k,s)*r(i)*r(j) + 2.d0 * v(k)
     $                    * (r(i)*djs + dis*r(j) )+ 2.d0 * v(s) * (r(i)
     $                    *djk + r(j)*dik)+ rrq * (dis*djk + djs*dik)

                     dy2s = dt2s*r(k) + t2*dks + ddt3ks

                     
                     ddt3ls = 2.d0 * (q(l,s)*r(i)*r(j) + v(l)*(dis*r(j)
     $                    + r(i)*djs) )
     $                    + 2.d0 * v(s) * (r(i)*djl + r(j)*dil)
     $                    + rrq * (dis*djl + djs*dil)

                     dt3s = 2.d0 * v(s)*r(i)*r(j) + rrq * (r(i)*djs +
     $                    r(j)*dis)
                     ddy3ls = ddt3ls*r(k) + dt3l*dks + dt3s * dkl

                     dz3s = dy2s*r(l) + y2*dls + ddy3ls
                     w2 = z2*r(s) + dz3s

c------------
c                     dt3s =
                     dy3s = dt3s*r(k) + t3 * dks

                     dz4s = dy3s*r(l) + y3*dls
                     w3 = z3*r(s) + dz4s

c-------------
                     w4 = z4*r(s) 
                     
                     d5q(i,j,k,l,s) = w1*r9 + w2*r11 + w3*r13 + w4
     $                    *r15
                     
                  end do
               end do
            end do
         end do
      end do


      return
      end
