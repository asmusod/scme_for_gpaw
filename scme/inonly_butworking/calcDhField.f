c----------------------------------------------------------------------+
c     Derivatives of the hexadecapole field                            |
c----------------------------------------------------------------------+
      subroutine dHpole(h, r, d1h, d2h, d3h, d4h, d5h)

      implicit none
      real*8 h(3,3,3,3), r(3), d1h(3), d2h(3,3), d3h(3,3,3)
      real*8 d4h(3,3,3,3), d5h(3,3,3,3,3) 

      integer delta, dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
      integer i, j, k, l, s
      real*8 r2, r3, r5, r7, r9, r11, r13, r15, r17, r19
      real*8 t1, t2, t3, y1, y2, y3, y4, z1, z2, z3, z4, z5, w1, w2, w3,
     $     w4, w5 
      real*8 r4h, v(3), g(3,3), d(3,3,3)
      
      real*8 dt1k, dt2k, dt3k
      real*8 dy1l, dt1l, ddt2kl, dy2l, dt2l, ddt3kl, dy3l, dt3l, dy4l
      real*8 ddt1ls, dt1s, dddt2kls, ddy2ls, dy1s, dz2s, ddt2ks, dy2s,
     $     dt2s, ddt2ls, dddt3kls, ddy3ls, dz3s, ddt3ks, dy3s, dt3s,
     $     ddt3ls, ddy4ls, dz4s, dy4s, dz5s
      

      r2 = r(1)**2 + r(2)**2 + r(3)**2
      r3 = sqrt(r2) * r2
      r5 = r3 * r2
      r7 = r5 * r2
      r9 = r7 * r2
      r11 = r9 * r2
      r13 = r11 * r2
      r15 = r13 * r2
      r17 = r15 * r2
      r19 = r17 * r2

      r9  =        1.d0 / r9
      r11 =       -9.d0 / r11
      r13 =       99.d0 / r13
      r15 =    -1287.d0 / r15
      r17 =    19305.d0 / r17
      r19 =  -328185.d0 / r19

      r4h = 0.d0
      do l = 1, 3
         v(l) = 0.d0
         do k = 1, 3
            g(l,k) = 0.d0
            do j = 1, 3
               d(l,k,j) = 0.d0
               do i = 1, 3
                  r4h    = r4h + h(i,j,k,l) * r(i) * r(j) * r(k) * r(l)
                  v(l)   = v(l) + h(i,j,k,l) * r(i) * r(j) * r(k) 
                  g(l,k) = g(l,k) + h(i,j,k,l) * r(i) * r(j)
                  d(l,k,j) = d(l,k,j) + h(i,j,k,l) * r(i)
               end do
            end do
         end do
      end do


      do i = 1, 3
         d1h(i) = (4.d0 * v(i)) * r9 + (r4h * r(i)) * r11

         do j = i, 3
            dij = delta(i,j)

c     2nd derivative
            t1 = 12.d0 * g(i,j)
            t2 = 4.d0 * (v(i)*r(j) + v(j)*r(i)) + r4h * dij
            t3 = r4h * r(i)*r(j)

            d2h(i,j) = t1 * r9 + t2 * r11 + t3 * r13

            do k = j, 3
               dik = delta(i,k)
               djk = delta(j,k)
               
c     3rd derivative
               dt1k = 24.d0 * d(i,j,k)
               dt2k = 4.d0 * (v(i)*djk + v(j)*dik + v(k) * dij) + 12.d0
     $              * (g(i,k)*r(j) + g(j,k)*r(i))
               dt3k = 4.d0 * v(k)*r(i)*r(j) + r4h * (r(i)*djk + r(j)*dik
     $              )
               y1 = dt1k 
               y2 = t1*r(k) + dt2k
               y3 = t2*r(k) + dt3k
               y4 = t3*r(k)
               
               d3h(i,j,k) = y1*r9 + y2*r11 + y3*r13 + y4*r15

               do l = k, 3
                  dil = delta(i,l)
                  djl = delta(j,l)
                  dkl = delta(k,l)
               
c     4th derivative
                  dy1l = 24.d0 * h(i,j,k,l)
                  z1 = dy1l
c-------
                  dt1l = 24.d0 * d(i,j,l)
                  ddt2kl = 12.d0 * (g(i,l)*djk + g(j,l)*dik + g(k,l)
     $                 * dij) + 24.d0 * (d(i,k,l)*r(j) + d(j,k,l)*r(i))
     $                 + 12.d0 * (g(i,k)*djl + g(j,k)*dil)
                  dy2l = dt1l*r(k) + t1*dkl + ddt2kl
                  z2 = y1*r(l) + dy2l

c-------
                  dt2l = 4.d0 * (v(i)*djl + v(j)*dil + v(l) * dij) + 12
     $                 .d0* (g(i,l)*r(j) + g(j,l)*r(i))
               
                  ddt3kl = 4.d0 * (3.d0 * g(k,l)*r(i)*r(j) 
     $                 + v(k) * (dil*r(j) + r(i)*djl) ) 
     $                 + r4h * (dil*djk + djl*dik) 
     $                 + 4.d0 * v(l) * (r(i)*djk + r(j)*dik)
                  dy3l = dt2l*r(k) + t2*dkl + ddt3kl
                  z3 = y2*r(l) + dy3l
                  
c-------
                  dt3l = 4.d0 * v(l)*r(i)*r(j) + r4h * (r(i)*djl + r(j)
     $                 *dil)
                  dy4l = dt3l*r(k) + t3 * dkl
                  z4 = y3*r(l) + dy4l

                  z5 = y4*r(l)
                  d4h(i,j,k,l) =  z1*r9 + z2*r11 + z3*r13 + z4*r15 +
     $                 z5*r17 
                  
                  do s = l, 3
                     dis = delta(i,s)
                     djs = delta(j,s)
                     dks = delta(k,s)
                     dls = delta(l,s)

c     5th derivative
                     ddt1ls = 24.d0 * h(i,j,l,s)
                     dt1s = 24.d0 * d(i,j,s)

                     dddt2kls = 24.d0 * (d(i,l,s)*djk + d(j,l,s)*dik +
     $                    d(k,l,s) * dij) + 24.d0 * (h(i,k,l,s)*r(j) +
     $                    h(j,k,l,s)*r(i) + d(i,k,l)*djs + d(j,k,l)*dis)
     $                    + 24.d0 * (d(i,k,s) *djl + d(j,k,s)*dil) 

                     
                     ddy2ls = ddt1ls*r(k) + dt1l*dks + dt1s*dkl +
     $                    dddt2kls

                     dy1s = 24.d0 * h(i,j,k,s)
                     dz2s = y1*dls + dy1s*r(l) + ddy2ls
                     w1 = z1*r(s) + dz2s
c----------

                     ddt2ks = 12.d0 * (g(i,s)*djk + g(j,s)*dik + g(k,s)
     $                    * dij) + 24.d0 * (d(i,k,s)*r(j) + d(j,k,s)*r(i
     $                    ))+ 12.d0 * (g(i,k)*djs + g(j,k)*dis)
                     dy2s = dt1s*r(k) + t1*dks + ddt2ks

                     dt2s = 4.d0 * (v(i)*djs + v(j)*dis + v(s) * dij) +
     $                    12.d0* (g(i,s)*r(j) + g(j,s)*r(i))

                     ddt2ls = 12.d0 * (g(i,s)*djl + g(j,s)*dil + g(l,s)
     $                    * dij) + 24.d0 * (d(i,l,s)*r(j) + d(j,l,s)*r(i
     $                    ))+ 12.d0 * (g(i,l)*djs + g(j,l)*dis)

                     dddt3kls = 12.d0 * (2.d0 * d(k,l,s)*r(i)*r(j) + g(k
     $                    ,l) * (dis*r(j) + r(i)*djs) )
     $                    + 4.d0 * (3.d0 * g(k,s) * (dil*r(j) + r(i)*djl
     $                    )+ v(k) * (dil*djs + dis*djl) )+ 4.d0 * v(s)
     $                    * (dil*djk + djl*dik)+ 4.d0 * (v(l) * (dis*djk
     $                    + djs*dik) + 3.d0* g(l,s) * (r(i)*djk + r(j)
     $                    *dik))


                     ddy3ls = ddt2ls*r(k) + dt2l*dks + dt2s*dkl +
     $                    dddt3kls

                     dz3s = dy2s*r(l) + y2*dls + ddy3ls

                     w2 = z2*r(s) + dz3s

c---------
c                    dt2s =

                     ddt3ks = 4.d0 * (3.d0 * g(k,s)*r(i)*r(j) 
     $                    + v(k) * (dis*r(j) + r(i)*djs) ) 
     $                    + r4h * (dis*djk + djs*dik) 
     $                    + 4.d0 * v(s) * (r(i)*djk + r(j)*dik)

                     dy3s = dt2s*r(k) + t2*dks + ddt3ks

                     dt3s = 4.d0 * v(s)*r(i)*r(j) + r4h * (r(i)*djs +
     $                    r(j)*dis)

                     ddt3ls = 4.d0 * (3.d0 * g(l,s)*r(i)*r(j) 
     $                    + v(l) * (dis*r(j) + r(i)*djs) ) 
     $                    + r4h * (dis*djl + djs*dil) 
     $                    + 4.d0 * v(s) * (r(i)*djl + r(j)*dil)

                     ddy4ls = ddt3ls*r(k) + dt3l*dks + dt3s * dkl

                     dz4s = dy3s*r(l) + y3*dls + ddy4ls

                     w3 = z3*r(s) + dz4s

c-------------
                     dy4s = dt3s*r(k) + t3 * dks
                     dz5s = dy4s*r(l) + y4*dls

                     w4 = z4*r(s) + dz5s

                     w5 = z5*r(s)

                     d5h(i,j,k,l,s) = w1*r11 + w2*r13 + w3*r15 + w4
     $                    *r17+w5*r19
                     
                  end do
               end do
            end do
         end do
      end do


      return
      end
