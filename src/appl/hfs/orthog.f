*
*     ------------------------------------------------------------------
*       O R T H O G
*     ------------------------------------------------------------------
*
      SUBROUTINE ORTHOG(LET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/SIGNF /SIGNFA
*
*     THIS SUBROUTINE CHECKS FOR POSSIBLE ORTHOGONALITY DUE TO
*     COUPLING DIFFERENCES OR UNEVEN PARITY
*
  101 FORMAT(37H DIFFERING RESULTANT ANGULAR MOMENTUM)
  102 FORMAT(52H ORTHOGONALITY IN COUPLING SCHEMES OF CONFIGURATIONS)
  103 FORMAT(59H THE TWO CONFIGURATIONS HAVE DIFFERING NUMBERS OF ELECTR
     :ONS)
  104 FORMAT(51H THE TWO CONFIGURATIONS HAVE DIFFERING TOTAL PARITY)
*
* --- DO PSI AND PSIP CONTAIN THE SAME NUMBERS OF ELECTRONS
*     DO PSI AND PSIP HAVE THE SAME TOTAL PARITY
*
      N5=0
      N6=0
      N7=0
      DO 89 I=1,IHSH
      L1=LJ(I)
      L2=NOSH1(I)
      L3=NOSH2(I)
      N5=N5+L2
      N6=N6+L3
      N7=N7+L1*(L2-L3)
   89 CONTINUE
*
*     CHECK ON NUMBER OF ELECTRONS
*
      IF (N5-N6) 21,22,21
   21 CONTINUE
   28 WRITE(IWRITE,103)
      GO TO 11
*
*     CHECK ON PARITY
*
   22 IF(N7-N7/2*2) 23,24,23
   23 CONTINUE
   25 WRITE(IWRITE,104)
      GO TO 11
   24 N72=N7/2
      SIGNFA=1.D0
      IF( (N72-(N72/2)*2).NE.0 ) SIGNFA=-1.D0
      N1=2*IHSH-1
      N2=IHSH+1
      N3=IHSH-1
      N4=IHSH-2
*
* --- IS THE FINAL STATE THE SAME FOR PSI AND PSIP
*
      DO 1 K=2,3
      IF(J1QN1(N1,K)-J1QN2(N1,K))2,1,2
    1 CONTINUE
      GO TO 3
    2 CONTINUE
   26 WRITE(IWRITE,101)
*
* --- THE TWO CONFIGURATIONS WILL HAVE ZERO HAMILTONIAN MATRIX ELEMENT
*
   11 LET=0
      RETURN
*
* --- COUPLING ORTHOGONALITY TEST FOR FIRST TWO SHELLS
*
    3 CONTINUE
   12 LET=1
      RETURN
      END