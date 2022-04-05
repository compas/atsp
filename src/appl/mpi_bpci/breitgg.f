*
*     ------------------------------------------------------------------
*       B R E I T G G
*     ------------------------------------------------------------------
*
      SUBROUTINE BREITGG(NEW,NZERO,IFIRST,idg,irfst,skip,nze,plnze)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irfst(*)
      PARAMETER (NTERMD=31)
      PARAMETER (MXIHSH=16)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3),iflag(3)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      common/diAGNL/IDIAG,JA,JB
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     ;             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),
     :      J1QN1(31,3),J1QN2(31,3),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QACMULT,ACMULT(1))
      COMMON /SPORB/ QACMULT
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/IQLSP,index(ntermd),IJK,termsh(ntermd),nterm
      pointer(pind,ind(nt)),(pnze_term,nze_term(njv)),
     :       (pnze_col,nze_col(njv)),(pnze_t,nze_t(njv)),
     :       (pnze_c,nze_c(njv))
      common/im_nze/pflsj,pnze_term,pnze_col,pind,nt,mxj,
     :              mnj,njv
      common/hmx_nze/pnze_t,pnze_c
!!!      common/NZZEL/nzzeall,nzznzero
      pointer(pflsj,flsj(nt,nt,2,njv))
!      integer nze_term(njv),nze_col(njv),ind(nt),col_tmp(njv)
      integer col_tmp(njv)
      integer nterm,mxj,mnj,njv,index
      double precision flsj
      COMMON/STEGG/IX,IGGG,IRHO,ISIG,IRHOP,ISIGP
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON /BREIT/ ISPORB,ISOORB,ISPSPN

      LOGICAL INCL,skip
      INTEGER ipos(3)
      pointer (plnze,lnze(nt,njv))
      nterm = nt;
*
* --- THIS PROGRAMME CONSIDERS EITHER SUPERPOSITION OF CONFIGURATIONS OR
*     MULTI-CONFIGURATIONAL HARTREE-FOCK WAVE FUNCTIONS.  USING THE
*     RESULT THAT THE TWO-ELECTRON HAMILTONIAN MATRIX ELEMENT
*     (PSI/V/PSIP)  CAN BE WRITTEN AS A SUM OF SLATER INTEGRALS, THE
*     PRESENT CODE  -  WEIGHTS  -  CALCULATES THE COEFFICIENTS OF THESE
*     INTEGRALS.  PSI AND PSIP ARE ALLOWED TO RUN OVER NCFG CONFIGURATNS
*
*
* --- CONSIDER (PSI/V/PSIP) AS PSI AND PSIP RUN OVER ALL CONFIGURATIONS
*

!      NFIRST = NCFG - NEW + 1
!      irow = max(jb,nfirst)
!      if (jb .gt. nzero) then
!	last = jb
!      else
!	last = ncfg
!      end if
!
      NFIRST = NCFG - NEW + 1
      irow = max(jb,nfirst)
      if (jb .gt. nzero) then
!       turn off the two-body relativistic operators (CFF: March_14,2001)
        isoorb = 0
        ispspn = 0
        iorborb = 0
      end if
      last = ncfg

      m = last-irow+1
      do 10 i = 1,m
	h(i,1) = 0.d0
  10  continue
      if (irel .ne. 0) then
	do i = 1,m
 	  h(i,2) = 0.d0
	  h(i,3) = 0.d0
        end do
      end if
      nrow(1) = 0
      nrow(2) = 0
      nrow(3) = 0

*<<<<<<<>>>>>>>>
*    do all in the current column
      do ja = irow, last
	 do iel = 1,maxorb
	   acmult(iel) = 0.d0
         end do
         INCL = .TRUE.
*     IF (JB.LE.NZERO.AND.IRFST(JB).EQ.1) INCL = .TRUE.
         IFLAG(1) = 0
         IFLAG(2) = 0
         IFLAG(3) = 0
         IDIAG = 0
         IF(JA.EQ.JB) THEN
            IDIAG=1
            IF(IDG.EQ.1) INCL = .TRUE.
         ENDIF
*
* ... Set up defining quantum numbers for each matrix element.
*
         IF(INCL) THEN
            CALL SHELLS(JA,JB,LET)
            IF(LET.NE.0) THEN
                IF(IHSH.GT.MXIHSH) STOP
*
* ... TEST ON POSSIBLE RECOUPLING ORTHOGONALITY.
*
                CALL SETUPGG
                IF(IX.LE.4) THEN
                   CALL ORTHOG(LET,INCL)
                   IF(LET.NE.0) THEN
                     IF (IFULL.NE.0) WRITE(IWRITE,77)
                     IF(IBUG1.NE.0 .OR. IBUG2.NE.0) CALL VIJOUT(JA,JB)
                     CALL NONBP
                   ENDIF
                ENDIF
            ENDIF
         ENDIF
      end do 
*      end column 

      m1 = nrow(1);
      m2 = nrow(2);
      m3 = nrow(3);

      nze = nze + m1
      write(iouhn) jb,m1,(h(i,1),i=1,m1),(jan(i,1),i=1,m1)
      if (.not. skip) then
        m2 = nrow(2)
        write(iouhz) jb,m2,(h(i,2),i=1,m2),(jan(i,2),i=1,m2)
        write(iouhs) jb,m3,(h(i,3),i=1,m3),(jan(i,3),i=1,m3)
      end if 

!   find which term is jb
      jb_term = 1; 
      do while (ind(jb_term) < jb)
         jb_term = jb_term + 1; 
!         if (jb_term == nt) exit;
      end do

      col_tmp(1:njv) = 0
***
*<<<<<<<<>>>>>>>>> find non zero
      jc = 1;
      do jj = mxj,mnj,-2
         do ic1 = 1, nt
           call lsp_cv(lsp(ic1),lli,ksi)
           if(jj>=abs(lli-ksi).and.jj<=(lli+ksi)) then
              lsp(ic1) = 2*(lsp(ic1)/2) + 1;
           else
              lsp(ic1) = 2*(lsp(ic1)/2);
           end if;
         end do
         if (mod(lsp(jb_term),2) == 0) then
            col_tmp(jc) = col_tmp(jc) + 1 
            nze_term(jc) = nze_term(jc) + 1
         else 
            ipos(1:3) = 0; 
            call advance(1,m1,ipos,ncfg)
            call advance(2,m2,ipos,ncfg)
            call advance(3,m3,ipos,ncfg)
            iipos = min(ipos(1),ipos(2),ipos(3))
            do while (iipos.le.ncfg)
               col_tmp(jc) = col_tmp(jc) + 1
               nze_term(jc) = nze_term(jc) + 1
               irow = min(nrow(1),nrow(2),nrow(3))
               if (nrow(1) .eq. irow) call advance(1,m1,ipos,ncfg)
               if (nrow(2) .eq. irow) call advance(2,m2,ipos,ncfg)
               if (nrow(3) .eq. irow) call advance(3,m3,ipos,ncfg)
               iipos = min(ipos(1),ipos(2),ipos(3)) 
            end do
         end if
         jc = jc + 1;
       end do
!       print*, nze_term(1), nze_term(2), nze_term(3)
       
       do jc = 1, njv
          if (col_tmp(jc) > nze_col(jc) ) nze_col(jc) = col_tmp(jc) 
       end do
*<<<<>>>>
**
!      write(0,'(5F14.8)') (h(i,1),i=1,m1)
!      write(0,'(5i14)') (jan(i,1),i=1,m1)
!      if (.not. skip) then
!      write(0,'(5F14.8)') (h(i,2),i=1,m2)
!      write(0,'(5i14)') (jan(i,2),i=1,m2)
!      write(0,'(5F14.8)') (h(i,3),i=1,m3)
!      write(0,'(5i14)') (jan(i,3),i=1,m3)
!      end if
*
   77 FORMAT(///30X,'MULTIPLYING FACTOR',11X,'TYPE OF INTEGRAL')
      END
