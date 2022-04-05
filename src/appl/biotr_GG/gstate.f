*
*     ------------------------------------------------------------------
*	G S T A T E
*     ------------------------------------------------------------------
*
      SUBROUTINE gstate(NFIRST,NLAST)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128)
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      CHARACTER*1 JCQN(15),JQNG(15),J3QNG(15)
      DIMENSION J1QN(15),J2QN(15),J3QN(15)
      CHARACTER LABEL(2)*8, LINE*72
      DATA LABEL/'INITIAL','FINAL'/
*
*      DATA DEFINING THE STATE IS READ IN AND PRINTED OUT.
*
    5 FORMAT(8(1X,A3,1H(,I2,1H)))
   55 FORMAT(8(1X,A3,1X,I2,1X))
CGG    6 FORMAT(15(1X,I1,A1,I1))
   36 FORMAT(15(1X,A1,A1,A1))
   24 FORMAT(//31H INITIAL STATE CONFIGURATIONS:-)
   25 FORMAT(/5H     ,I5,3H.  ,10(1X,A3,1H(,I2,1H)))
   26 FORMAT(11X,10(1X,4X,I1,A1,I1))
CGG   27 FORMAT(22X,9(1X,4X,I1,A1,I1))
   27 FORMAT(22X,9(1X,4X,A1,A1,I1))
   28 FORMAT(  31H ----------------------------  /)
   29 FORMAT(//29H FINAL STATE CONFIGURATIONS:-)
   30 FORMAT(2X,'ELECTRON ',A3,' NOT FOUND IN THE LIST OF ELECTRONS',
     :   ' FOR THE ',A8,' STATE')
      IF (NFIRST .EQ. 1) THEN
C        WRITE(IWRITE,24)
         IRD = iuc(1)
        ELSE
C        WRITE(IWRITE,29)
         IRD = iuc(2)
      END IF
      WRITE(IWRITE,28)
      DO 2 I=NFIRST,NLAST
         READ(IRD,'(A72)') LINE
         N=0
         J=2
   60    IF (LINE(J:J+2) .NE. '   ' .and. N .LT. (8)) THEN
	 N = N + 1
	 J = J +8
	 GO TO 60
      END IF
      NOCCSH(I) = N
      READ(LINE,55)       (NOCORB(J,I),NELCSH(J,I),J=1,N)
      K=I
      IF(NFIRST.NE.1) K=I-NFIRST+1
Cdbg  WRITE(IWRITE,25) K,(NOCORB(J,I),NELCSH(J,I),J=1,N)
      NCOM1 = NCOM + 1
      NOR11 = NCOM1 + NORBI
      DO 61 J=1,N
      DO 63 JJ = 1,MAXORB
      IF (NFIRST .EQ. 1 .AND. JJ .GE. NOR11) GO TO 65
      IF(NFIRST .NE. 1 .AND. JJ .GE. NCOM1 .AND. JJ .LT. NOR11) GO TO 63
      IF(NOCORB(J,I).EQ.IAJCMP(JJ)) THEN
         NOCORB(J,I) = JJ
         GO TO 61
      END IF
   63    CONTINUE
*
*        ELECTRON NOT FOUND IN THE LIST
*
   65 WRITE(IWRITE,30) NOCORB(J,I),LABEL(NFIRST)
      STOP
   61 CONTINUE
      M=2*N-1
      N1=N+1
CGG      READ(IRD,6)      (J3QN(J),JCQN(J),J1QN(J),J=1,M)
      READ(IRD,36)     (J3QNG(J),JCQN(J),JQNG(J),J=1,M)
Cdbg  WRITE(IWRITE,26) (J3QNG(J),JCQN(J),J1QN(J),J=1,N)
Cdbg      IF(N.EQ.1) GO TO 64
Cdbg  WRITE(IWRITE,27) (J3QN(J),JCQN(J),J1QN(J),J=N1,M)
Cdbg   64 CONTINUE
      DO 62 J=1,M
      J2QN(J) = 2*LVAL(JCQN(J)) + 1
CGG
      J1QN(J) = NUMVAL(JQNG(J))
      J3QN(J) = ICHAR(J3QNG(J)) - ICHAR('1') + 1
CGG
   62 J1QNRD(J,I)= (J3QN(J)*64 + J2QN(J))*64 + J1QN(J)
    2 CONTINUE
      RETURN
      END
