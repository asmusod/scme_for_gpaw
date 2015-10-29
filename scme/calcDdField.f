c----------------------------------------------------------------------+
c     Derivatives of the dipole field                                  |
c----------------------------------------------------------------------+
      subroutine dDpole(d, r, d1d, d2d, d3d, d4d, d5d)

      implicit none
      real*8 d(3), r(3), d1d(3), d2d(3,3), d3d(3,3,3)
      real*8 d4d(3,3,3,3), d5d(3,3,3,3,3) 

      integer delta, dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
      integer i, j, k, l, s
      real*8 r2, r3, r5, r7, r9, r11, r13, rd
      real*8 t1, t2, y1, y2, y3, z1, z2, z3, w1, w2, w3, w4

      real*8 dt1k, dt2k, dt1l, ddt2kl, dy2l, dt2l, dy3l, dddt2kls
      real*8 dt1s, ddy2ls, dz1s, ddt2ks, dy2s, ddt2ls, dt2s
      real*8 ddy3ls, dz2s, dy3s, dz3s

      r2 = r(1)**2 + r(2)**2 + r(3)**2
      r3 = sqrt(r2) * r2
      r5 = r3 * r2
      r7 = r5 * r2
      r9 = r7 * r2
      r11 = r9 * r2
      r13 = r11 * r2

      r3  =        1.d0 / r3
      r5  =       -3.d0 / r5
      r7  =       15.d0 / r7
      r9  =     -105.d0 / r9
      r11 =      945.d0 / r11
      r13 =   -10395.d0 / r13

      rd = r(1)*d(1) + r(2)*d(2) + r(3)*d(3)

      do i = 1, 3
         d1d(i) = d(i) * r3 + rd * r(i) * r5

         do j = i, 3
            dij = delta(i,j)

c     2nd derivative
            t1 = d(i)*r(j) + d(j)*r(i) + rd * dij
            t2 = rd * r(i)*r(j)

            d2d(i,j) = t1 * r5 + t2 * r7 

            do k = j, 3
               dik = delta(i,k)
               djk = delta(j,k)
               
c     3rd derivative

               dt1k = d(i)*djk + d(j)*dik + d(k) * dij
               y1 = dt1k
               dt2k = d(k) * r(i)*r(j) + rd * (r(i)*djk + dik*r(j))
               y2 = t1*r(k) + dt2k
               y3 = t2 * r(k)
               
               d3d(i,j,k) = y1*r5 + y2*r7 + y3*r9 

               do l = k, 3
                  dil = delta(i,l)
                  djl = delta(j,l)
                  dkl = delta(k,l)
               
c     4th derivative

                  dt1l = d(i)*djl + d(j)*dil + d(l) * dij
                  ddt2kl = d(k) * (r(i)*djl + dil*r(j))
     $                 + d(l) * (r(i)*djk + dik*r(j))
     $                 + rd * (dil*djk + dik*djl)    
                  dy2l = dt1l*r(k) + t1*dkl + ddt2kl
                  z1 = y1*r(l) + dy2l

                  dt2l = d(l) * r(i)*r(j) + rd * (r(i)*djl + dil*r(j))
                  dy3l = dt2l * r(k) + t2 * dkl
                  z2 = y2*r(l) + dy3l

                  z3 = y3*r(l)

                  d4d(i,j,k,l) =  z1*r7 + z2*r9 + z3*r11
                  
                  do s = l, 3
                     dis = delta(i,s)
                     djs = delta(j,s)
                     dks = delta(k,s)
                     dls = delta(l,s)

c     5th derivative
                     dddt2kls = d(k) * (dis*djl + dil*djs)
     $                    + d(l) * (dis*djk + dik*djs)
     $                    + d(s) * (dil*djk + dik*djl)    
                     dt1s = d(i)*djs + d(j)*dis + d(s) * dij

                     ddy2ls = dt1l*dks + dt1s*dkl + dddt2kls
                     dz1s = y1*dls + ddy2ls
                     w1 = dz1s


                     ddt2ks = d(k) * (r(i)*djs + dis*r(j))
     $                    + d(s) * (r(i)*djk + dik*r(j))
     $                    + rd * (dis*djk + dik*djs)    
                     dy2s = dt1s*r(k) + t1*dks + ddt2ks

                     ddt2ls = d(l) * (dis*r(j) + r(i)*djs)
     $                    + d(s) * (r(i)*djl + dil*r(j))
     $                    + rd * (dis*djl + dil*djs)
                     dt2s = d(s) * r(i)*r(j) + rd * (r(i)*djs + dis*r(j
     $                    ))
                     ddy3ls = ddt2ls * r(k) + dt2l * dks + dt2s * dkl
                     
                     dz2s = dy2s*r(l) + y2*dls + ddy3ls
                     w2 = z1 * r(s) + dz2s

                     dy3s = dt2s * r(k) + t2 * dks
                     dz3s = dy3s*r(l) + y3*dls
                     w3 = z2 * r(s) + dz3s

                     w4 = z3 * r(s)
                     
                     d5d(i,j,k,l,s) = w1*r7 + w2*r9 + w3*r11 + w4
     $                    *r13
                     
                  end do
               end do
            end do
         end do
      end do


      return
      end
