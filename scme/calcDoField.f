c----------------------------------------------------------------------+
c     Derivatives of the octopole field                                |
c----------------------------------------------------------------------+
      subroutine dOpole(o, r, d1o, d2o, d3o, d4o, d5o)

      implicit none
      real*8 o(3,3,3), r(3), d1o(3), d2o(3,3), d3o(3,3,3)
      real*8 d4o(3,3,3,3), d5o(3,3,3,3,3) 

      integer delta, dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
      integer i, j, k, l, s
      real*8 r2, r3, r5, r7, r9, r11, r13, r15, r17
      real*8 t1, t2, t3, y1, y2, y3, y4, z1, z2, z3, z4, w1, w2, w3, w4,
     $     w5 
      real*8 r3o, v(3), g(3,3)

      real*8 dt1k, dt2k, dt3k
      real*8 dt1l, ddt2kl, dy2l, dt2l, ddt3kl, dy3l, dt3l, dy4l
      real*8 dt1s, dddt2kls, ddy2ls, dz1s, ddt2ks, dy2s, ddt2ls, dt2s,
     $     dddt3kls, ddy3ls, dz2s, ddt3ks, dy3s, dt3s, ddt3ls, ddy4ls,
     $     dz3s, dy4s, dz4s


      r2 = r(1)**2 + r(2)**2 + r(3)**2
      r3 = sqrt(r2) * r2
      r5 = r3 * r2
      r7 = r5 * r2
      r9 = r7 * r2
      r11 = r9 * r2
      r13 = r11 * r2
      r15 = r13 * r2
      r17 = r15 * r2

      r7  =       1.d0 / r7
      r9  =      -7.d0 / r9
      r11 =      63.d0 / r11
      r13 =    -693.d0 / r13
      r15 =    9009.d0 / r15
      r17 = -135135.d0 / r17

      r3o = 0.d0
      do k = 1, 3
         v(k) = 0.d0
         do j = 1, 3
            g(j,k) = 0.d0
            do i = 1, 3
               r3o    = r3o + o(i,j,k) * r(i) * r(j) * r(k)
               v(k)   = v(k) + o(i,j,k) * r(i) * r(j) 
               g(j,k) = g(j,k) + o(i,j,k) * r(i)
            end do
         end do
      end do

      do i = 1, 3
         d1o(i) = (3.d0 * v(i)) * r7 + (r3o * r(i)) * r9

         do j = i, 3
            dij = delta(i,j)

c     2nd derivative
            t1 = 6.d0 * g(i,j)
            t2 = 3.d0 * (v(i)*r(j) + v(j)*r(i)) + r3o * dij
            t3 = r3o * r(i)*r(j)

            d2o(i,j) = t1 * r7 + t2 * r9 + t3 * r11

            do k = j, 3
               dik = delta(i,k)
               djk = delta(j,k)
               
c     3rd derivative
               dt1k = 6.d0 * o(i,j,k)
               dt2k = 3.d0 * (v(i)*djk + v(j)*dik + v(k) * dij) + 6.d0
     $              * (g(i,k)*r(j) + g(j,k)*r(i))
               dt3k = 3.d0 * v(k)*r(i)*r(j) + r3o * (r(i)*djk + r(j)*dik
     $              )
               y1 = dt1k 
               y2 = t1*r(k) + dt2k
               y3 = t2*r(k) + dt3k
               y4 = t3*r(k)
               
               d3o(i,j,k) = y1*r7 + y2*r9 + y3*r11 + y4*r13

               do l = k, 3
                  dil = delta(i,l)
                  djl = delta(j,l)
                  dkl = delta(k,l)
               
c     4th derivative
                  dt1l = 6.d0 * o(i,j,l)
                  ddt2kl = 6.d0 * (g(i,l)*djk + g(j,l)*dik + g(k,l)
     $                 * dij) + 6.d0 * (o(i,k,l)*r(j) + o(j,k,l)*r(i))
     $                 + 6.d0 * (g(i,k)*djl + g(j,k)*dil)
                  dy2l = dt1l*r(k) + t1*dkl + ddt2kl
                  z1 = y1*r(l) + dy2l

c-------
                  dt2l = 3.d0 * (v(i)*djl + v(j)*dil + v(l) * dij) + 6
     $                 .d0* (g(i,l)*r(j) + g(j,l)*r(i))
               
                  ddt3kl = 3.d0 * (2.d0 * g(k,l)*r(i)*r(j) 
     $                 + v(k) * (dil*r(j) + r(i)*djl) ) 
     $                 + r3o * (dil*djk + djl*dik) 
     $                 + 3.d0 * v(l) * (r(i)*djk + r(j)*dik)
                  dy3l = dt2l*r(k) + t2*dkl + ddt3kl
                  z2 = y2*r(l) + dy3l
                  
c-------
                  dt3l = 3.d0 * v(l)*r(i)*r(j) + r3o * (r(i)*djl + r(j)
     $                 *dil)
                  dy4l = dt3l*r(k) + t3 * dkl
                  z3 = y3*r(l) + dy4l

                  z4 = y4*r(l)
                  d4o(i,j,k,l) =  z1*r9 + z2*r11 + z3*r13 + z4*r15

                  
                  do s = l, 3
                     dis = delta(i,s)
                     djs = delta(j,s)
                     dks = delta(k,s)
                     dls = delta(l,s)


c     5th derivative
                     dt1s = 6.d0 * o(i,j,s)

                     dddt2kls = 6.d0 * (o(i,l,s)*djk + o(j,l,s)*dik +
     $                    o(k,l,s) * dij + o(i,k,l)*djs + o(j,k,l)*dis + 
     $                    o(i,k,s)*djl + o(j,k,s)*dil)

                     ddy2ls = dt1l*dks + dt1s*dkl + dddt2kls
           
                     dz1s = y1*dls + ddy2ls
                     w1 = dz1s
c----------

                     ddt2ks = 6.d0 * (g(i,s)*djk + g(j,s)*dik + g(k,s)
     $                    * dij) 
     $                    + 6.d0 * (o(i,k,s)*r(j) + o(j,k,s)*r(i))
     $                    + 6.d0 * (g(i,k)*djs + g(j,k)*dis)
                     dy2s = dt1s*r(k) + t1*dks + ddt2ks

                     ddt2ls = 6.d0 * (g(i,s)*djl + g(j,s)*dil + g(l,s)
     $                    * dij) 
     $                    + 6.d0 * (o(i,l,s)*r(j) + o(j,l,s)*r(i))
     $                    + 6.d0 * (g(i,l)*djs + g(j,l)*dis)

                     dt2s = 3.d0 * (v(i)*djs + v(j)*dis + v(s) * dij) +
     $                    6.d0* (g(i,s)*r(j) + g(j,s)*r(i))

                     dddt3kls = 6.d0 * (o(k,l,s)*r(i)*r(j) + g(k,l)
     $                    * (dis*r(j) + r(i)*djs))
     $                    + 6.d0 * g(k,s) * (dil*r(j) + r(i)*djl)
     $                    + 3.d0 * v(k) * (dil*djs + dis*djl)
     $                    + 3.d0 * v(s) * (dil*djk + djl*dik)
     $                    + 6.d0 * g(l,s) * (r(i)*djk + r(j)*dik)
     $                    + 3.d0 * v(l) * (dis*djk + djs*dik)

                     ddy3ls = ddt2ls*r(k) + dt2l*dks + dt2s*dkl +
     $                    dddt3kls

                     dz2s = dy2s*r(l) + y2*dls + ddy3ls

                     w2 = z1*r(s) + dz2s

c------------
c                     dt2s = 

                     ddt3ks = 6.d0 * g(k,s)*r(i)*r(j) + 3.d0 * v(k)
     $                    * (dis*r(j) + r(i)*djs)
     $                    + 3.d0 * v(s) * (r(i)*djk + r(j)*dik)
     $                    + r3o * (dis*djk + djs*dik)

                     dy3s = dt2s*r(k) + t2*dks + ddt3ks

                     dt3s = 3.d0 * v(s)*r(i)*r(j) + r3o * (r(i)*djs +
     $                    r(j)*dis)

                     ddt3ls = 6.d0 * g(l,s)*r(i)*r(j) + 3.d0 * v(l)
     $                    * (dis*r(j) + r(i)*djs)
     $                    + 3.d0 * v(s) * (r(i)*djl + r(j)*dil)
     $                    + r3o * (dis*djl + djs*dil)

                     ddy4ls = ddt3ls*r(k) + dt3l*dks + dt3s * dkl

                     dz3s = dy3s*r(l) + y3*dls + ddy4ls
                                       
                     w3 = z2*r(s) + dz3s

c-------------
                     dy4s = dt3s*r(k) + t3 * dks
                     dz4s = dy4s*r(l) + y4*dls

                     w4 = z3*r(s) + dz4s

c---------
                     w5 = z4*r(s)
                     
                     d5o(i,j,k,l,s) = w1*r9 + w2*r11 + w3*r13 + w4
     $                    *r15 + w5*r17
                     
                  end do
               end do
            end do
         end do
      end do


      return
      end
