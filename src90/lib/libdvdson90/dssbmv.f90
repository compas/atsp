!***********************************************************************
      SUBROUTINE DSSBMV(A, INDROW, INDCOL, IUPPER, NZER, TEMPB, TEMPC, N, M, B&
         , C) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!***********************************************************************
!
!     Computes the product of matrix A with a block of vectors B(N,M)
!               C=A B
!     where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
!     parts of the matrix are kept and this is managed by Indices
!
!     Subroutines called:
!     dgathr,dscatr,ddot,daxpy,dini
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
       
      USE dgathr_I 
      USE ddot_I 
      USE daxpy_I 
      USE dscatr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NZER 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      LOGICAL , INTENT(IN) :: IUPPER 
      INTEGER  :: INDROW(NZER) 
      INTEGER , INTENT(IN) :: INDCOL(N) 
      REAL(DOUBLE)  :: A(NZER) 
      REAL(DOUBLE)  :: TEMPB(N) 
      REAL(DOUBLE)  :: TEMPC(N) 
      REAL(DOUBLE)  :: B(N*M) 
      REAL(DOUBLE)  :: C(N*M) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOFFS, ISTART, ICOL, NUMELEM, ICUR, IV 
      REAL(DOUBLE) :: DIAG 
!-----------------------------------------------
 
!***********************************************************************
!
!   on entry
!   --------
!   N  the dimension (rank) of the matrix A
!   B      the matrix (block of vectors) to multiply A with
!   A  Linear array keeping the the nonzero elements of
!      the matrix A. It stores columns one after the other
!      starting from the first. It is either the upper
!      triangular part or the lower part depending on logical
!      iupper. (max elements that may contain= n^2/(2*nodes))
!   INDCOL It is the index showing for each column where that
!      column ends in array A.
!   INDROW It is the index showing for each element of A to its
!      row number in the square matrix
!   iupper logical. If .true. the upper part of A is used.
!   TEMPB,TEMPC Linear scratch arrays (dim=N)
!
!   on exit
!   -------
!
!   C   the result of the multiplication (dim=NxM)
!
!***********************************************************************
      IOFFS = 1 
      IF (IUPPER) IOFFS = 0 
 
      ISTART = 1 
      CALL DINIT (N*M, 0.D0, C, 1) 
 
      DO ICOL = 1, N 
         NUMELEM = INDCOL(ICOL) - ISTART + 1 
         ICUR = 1 
 
         DO IV = 1, M 
            CALL DGATHR (NUMELEM, B(ICUR), 1, INDROW(ISTART), 1, TEMPB, 1) 
            CALL DGATHR (NUMELEM, C(ICUR), 1, INDROW(ISTART), 1, TEMPC, 1) 
 
            DIAG = C(ICUR-1+ICOL) + DDOT(NUMELEM,A(ISTART),1,TEMPB,1) 
            CALL DAXPY (NUMELEM - 1, B(ICUR-1+ICOL), A(ISTART+IOFFS), 1, TEMPC(&
               IOFFS+1), 1) 
 
            CALL DSCATR (NUMELEM, TEMPC, 1, INDROW(ISTART), 1, C(ICUR), 1) 
            C(ICUR-1+ICOL) = DIAG 
            ICUR = ICUR + N 
         END DO 
 
         ISTART = INDCOL(ICOL) + 1 
      END DO 
 
      RETURN  
      END SUBROUTINE DSSBMV 
