c                                                                                         
c   Random number generator:                                                              

      FUNCTION RANb(IX,J)

      implicit real*8 (a-h,o-z)

      INTEGER*4 J,IX,I,I1,I2,I3,I4,I5,I6,I7,I8,I9
      DATA I2/16807/,I4/32768/,I5/65536/,I3/2147483647/

      I6=IX/I5
      I7=(IX-I6*I5)*I2
      I8=I7/I5
      I9=I6*I2+I8
      I1=I9/I4
      IX=(((I7-I8*I5)-I3)+(I9-I1*I4)*I5)+I1
      IF(IX.LT.0) IX=IX+I3
      RANb=(4.656613E-10)*FLOAT(IX)

      RETURN
      END
 
