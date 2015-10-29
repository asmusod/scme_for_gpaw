c   Dblexp:    June 92
c  This routine calculates FOFRHO
c  Generalized to arbitrary order

	SUBROUTINE FINT(rho1,natm1,ipot,emb)
	
	implicit real*8 (a-h,o-z)

        include '../commonblks/parameters.cmn'
        include '../commonblks/compotent.cmn'
	
        dimension emb(MAXATOMS,3),natm1(3),rho1(MAXATOMS,3)

c     ---------------------------------------------------------------

	do i = 1,natm1(ipot)
          emb(i,ipot) = 0.0
	end do
        

	RETURN
	END


