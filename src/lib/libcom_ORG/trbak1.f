*
*     --------------------------------------------------------------
*        T R B A K 1
*     --------------------------------------------------------------
*
*
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z) 
* 
      INTEGER I,J,K,L,M,N,NM 
      DOUBLE PRECISION A(NM,N),E(N),Z(NM,M) 
      DOUBLE PRECISION S 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1, 
*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). 
* 
*     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC 
*     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING 
*     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1. 
* 
*     ON INPUT: 
* 
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
*          DIMENSION STATEMENT; 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS- 
*          FORMATIONS USED IN THE REDUCTION BY  TRED1 
*          IN ITS STRICT LOWER TRIANGLE; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL 
*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY; 
* 
*        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED; 
* 
*        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED 
*          IN ITS FIRST M COLUMNS. 
* 
*     ON OUTPUT: 
* 
*        Z CONTAINS THE TRANSFORMED EIGENVECTORS 
*          IN ITS FIRST M COLUMNS. 
* 
*     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
      IF (M .EQ. 0) GO TO 200 
      IF (N .EQ. 1) GO TO 200 
* 
      DO 140 I = 2, N 
         L = I - 1 
         IF (E(I) .EQ. 0.0D0) GO TO 140 
* 
         DO 130 J = 1, M 
            S = 0.0D0 
* 
            DO 110 K = 1, L 
  110       S = S + A(I,K) * Z(K,J) 
*     :::::::::: DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1. 
*                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW :::::::::: 
            S = (S / A(I,L)) / E(I) 
* 
            DO 120 K = 1, L 
  120       Z(K,J) = Z(K,J) + S * A(I,K) 
* 
  130    CONTINUE 
* 
  140 CONTINUE 
* 
  200 RETURN 
*     :::::::::: LAST CARD OF TRBAK1 :::::::::: 
      END 
