************************************************************************
*
*
* --- This SUBROUTINE computes the Breit-Pauli matrix elements
*     for a variety of operators:
*
*     ONE-ELECTRON OPERATOR (L-integrals)
*     ELECTROSTATIC INTERACTION (F, G, R) integrals)
*     SPIN-ORBIT INTERACTION  (Z-integrals)
*     SPIN-OTHER-ORBIT INTERACTION (N, V integrals)
*     SPIN-SPIN INTERACTION (S-integrals)
*
*     The matrix is decomposed into three parts:
*       Hnr -- non-relativistic integrals 
*       HZ  -- contributions from Z, N, and V integrals
*       HS  -- contributions from S integrals.
*
*
*****:******************************************************************
*
      subroutine BREVAL(rl,skip,ec,nze,onlydvd,iounew,njv,maxj,minj,
     :                  PFLSJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NTERMD=31,LINT=500)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************
*
        logical onlydvd
      CHARACTER ANS*2, SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
      COMMON/PARAM/Hh,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON /BREIT/ ISPORB,ISOORB,ISPSPN
*****
*       this is for the in/out files to be known.
*       We have to keep different names for the routines below breval.
*       therefore we append an "x" at the end of each one that is used
*       with a different name.
      COMMON /INOUT/IREADx,IWRITEx,ISCWx,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
*****
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFGX
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QSIGNFA,SIGNFA(1))
      COMMON/PHASES/QSIGNFA,ICSTAS
      POINTER(QACMULT,ACMULT(1))
      COMMON /SPORB/ QACMULT
      LOGICAL rel,rl
      POINTER (qintptr,intptr(1)),(qpackn,ipackn(1)),(qvalue,value(1))
      COMMON /TABLE/qintptr,qpackn,qvalue,lmax
*     .. radial common
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/IQLSP,index(ntermd),IJK,termsh(ntermd),nterm
*
      POINTER(QIRFST,IRFST(1))
      DIMENSION NFLG(20)
      LOGICAL skip,lnze
      EXTERNAL INITT
      pointer(pflsj,nze_flsj(1)),(plnze,lnze(1,njv))
      double precision nze_flsj
*
*
  105 FORMAT (49H ISPORB=0 AND ISOORB=1 CAUSES THE PROGRAM TO FAIL,
     :  34H BECAUSE THE BLUME WATSON FORMULAE,/
     :  50H CANNOT BE USED FOR CLOSED SUBSHELLS.  TO OVERCOME,
     :  34H THIS, THE CODE HAS SET ISPORB = 1//)
   11 FORMAT(////' THE TYPE OF CALCULATION IS DEFINED BY ',
     :  'THE FOLLOWING PARAMETERS - '/
     : 5X,22H BREIT-PAULI OPERATORS,13X,8HIREL   =,I2/
     : 5X,27H PHASE CONVENTION PARAMETER,8X,8HICSTAS =,I2/)
   13 FORMAT(40H RELATIVISTIC OPERATORS TO BE INCLUDED -/5X,13H SPIN-ORB
     :IT (,I1,22H),  SPIN-OTHER-ORBIT (,I1,15H),  SPIN-SPIN (,I1,1H)/)
   24 FORMAT(36H0INITIAL DEBUG: IN 1-ELECTRON PART =,I2,2H ,,5X,
     : 20H IN 2-ELECTRON PART=,I2/16X,23HIN RECOUPLING PACKAGE =,I2/)
   42 FORMAT(//' MATRIX ELEMENTS CONSTRUCTED USING ',
     :       'THE SPHERICAL HARMONIC PHASE CONVENTION OF'/)
   43 FORMAT(/16X,47HCONDON AND SHORTLEY, THEORY OF ATOMIC STRUCTURE/
     :16X,47H-----------------------------------------------///)
   44 FORMAT(25X,42HFANO AND RACAH, IRREDUCIBLE TENSORIAL SETS/25X,42H--
     :----------------------------------------///)
   50 FORMAT(/20X,'======================='/
     :        20X,' B R E I T - P A U L I '/
     :        20X,'======================='/)
   78 FORMAT(19H DEBUG PARAMETERS -/5X,16H NBUG6(TENSOR) =,I2/5X,
     : ' NBUG7(RELATIVISTIC OPERATORS - SO,SOO,SS) =',I2//)
*
*     .. set unit numbers for "Breit" code (COMMON /INFORM/)
      iread = iuc
      iscw = ISCWx
      iwrite = IWRITEx
      iout =7
*     .. set rel in COMMON to be rl(input argument)
      rel = rl
      NIJ1 = 0
      NIJ2 = 0
      NIJ3 = 0
CG
      ISOTOP=0
      ICOLOM=0
      IORBORB=0
CG

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if (myid.eq.0) then
   10    WRITE(IWRITE,50)
      WRITE(ISCW,'(A/A/A/A)') ' Indicate the type of calculation ',
     : ' 0 => non-relativistic Hamiltonian only;',
     : ' 1 => one or more relativistic operators only;',
     : ' 2 => non-relativistic operators and selected relativistic:  '
         READ(5,*) IREL
         IF (IREL.NE.1) ICOLOM=1
         skip = .true.
         if (irel .ne. 0) then
           skip = .false.
         end if
         IFULL = 0
*     WRITE(ISCW,'(A)') ' Is full print-out requested? (Y/N) '
*     READ(5,'(A2)') ANS
*     IF (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') IFULL = 1
*
*     Determine basic parameters  (CS phases selected)
*
*     WRITE(ISCW,'(/A)')
*    :    ' Phases:- Condon and Shortley or Fano and Racah ? (CS/FR) '
*     READ(5,'(A2)') ANS
*        IF ( ANS .EQ. 'CS' .OR. ANS .EQ. 'cs') THEN
            ICSTAS = 1
*          ELSE
*           ICSTAS = 0
*        END IF
         ISPORB = 0
         ISOORB = 0
         ISPSPN = 0
         IORBORB=0
         IF (IREL .NE. 0) THEN
            ISPORB = 1
            ISOORB = 1
            ISPSPN = 1
            IORBORB=1
            WRITE(ISCW,*) ' All relativistic operators ? (Y/N) '
            READ(5,'(A2)') ANS
            IF ( ANS .EQ. 'N ' .OR. ANS .EQ. 'n ') THEN
               WRITE(ISCW,*) 'Spin-orbit,Spin-other-orbit',
     :         ',Spin-spin,Orbit-Orbit (0/1)'
               READ(5,*) ISPORB,ISOORB,ISPSPN,IORBORB
            END IF
            IF(ISPORB.EQ.0.AND.ISOORB.NE.0) THEN
              ISPORB = 1
              WRITE (IWRITE,105)
            END IF
         END IF
      endif
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            ICSTAS = 1
*
*  --:  Determine debug parameters
*
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
      NBUG6 = 0
      NBUG7 = 0
*
      WRITE(IWRITE,11) IREL,ICSTAS
*
*     SET FACTORIALS AND LOG OF FACTORIALS
*
       IF(IORBORB.EQ.1) MSOO=MSOO+3
       call MPI_BCAST(isporb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(isoorb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ispspn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(IORBORB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(irel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(skip,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(icolom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       ISPORBG = ISPORB
       ISOORBG = ISOORB
       ISPSPNG = ISPSPN
       IORBORBG = IORBORB
 79    CALL FACTRL(32)

*
* --- READ IN THE SET OF CONFIGURATIONS
*
*.      everybody reads the clist (cfgo1 etc)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

*	HOwever they read it in a sequential manner 
*	so that the disk is not thrashed by many accesses.
*

        CALL ACNFIG(NCLOSD);
        NCFG=NCFGX

      call alloc(qirfst,ncfg,4)
      IDG = 0
      NZERO = NCFG
      NEW = NCFG
      IRFST(1:ncfg) = 1
      ISTART = 1
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if (myid.eq.0) then
        WRITE(ISCW,'(A)') ' All Interactions? (Y/N): '
        READ (5,'(A2)') ANS
        IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
          WRITE(ISCW,'(A,I3,A/A)') ' Of the ',NCFG,
     :         ' configurations, how many at the end are new? ',
     :         ' How many configurations define the zero-order set?'
          READ (5,*) NEW,NZERO
          IF (NEW .EQ. 0) NEW = NCFG
          ISTART = NCFG - NEW + 1
          IF (NZERO .eq. 0) NZERO = NCFG
          IF (IREL .NE. 0 .AND. NZERO .NE. NCFG) THEN
             WRITE(ISCW,'(A)') 'Rel interaction with all the ',
     :                       'zero block?'
             READ (5,'(A2)') ANS
             IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
               DO 553 I = 1,(20)
553            NFLG(I) = 0
               WRITE(ISCW,*) ' Define your reference set : FORM(20I3)'
               READ (5,'(20I3)') (NFLG(I),I=1,(20))
               DO 551 I = 1,(ISTART-1)
               IRFST(I) = 0
               DO 552 J = 1,(20)
               IF (NFLG(J) .EQ. I) THEN
                  IRFST(I) = 1
                  GO TO 552
               ENDIF
552            CONTINUE
551            CONTINUE
            ENDIF
            WRITE(ISCW,*) ' Diagonal rel corrections ? (Y/N): '
            READ (5,'(A2)') ANS
            IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') IDG = 1
         ENDIF
         ISTRICT = 1
         WRITE(ISCW,*)  ' Restricted Two-body interactions? (Y/N); '
         READ (5,'(A2)') ANS
         IF (ANS .NE. 'Y' .and. ANS .NE. 'y' ) ISTRICT = 0
      END IF
      END IF
*       ..send imagnt and rest

*          ..sending /imagnt/ 1real = 2integers
!           len = 3*(2) + 8
!      call MPI_BCAST(const,len,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!      call MPI_BCAST(consts,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(new,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nzero,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nflg,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irfst,ncfg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icolom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


           
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* ...  Start the calculation.  Allocate memory for Tables of 
*      Integrals and for radial functions
*
      nwf = maxorb+nclosd
      call alctab(ncfg,nwf,skip)
*    
*     .. read radial functions, but first determint the number
*        of electrons
      nnel = 0
      do 65 i = 1,nclosd
        nnel = nnel + 2*(2*ljclsd(i)+1)
   65 continue
      do 66 i = 1,noccsh(1)
        nnel = nnel + nelcsh(i,1)
   66 continue
      call readw(rel,ec,nnel,skip,onlydvd,iounew)
*
* If only a davidson procedure is to be performed the creation 
* of breita is all skipped.
*
        if (onlydvd) goto 291

*     .. find maximum l-value
      lmax = 0
      do i=1,maxorb
         lmax = max0(lmax,ljcomp(i))
      end do
*     
*     .. generate list of configurations
      call genintbr(nclosd,maxorb,lmax,qljcomp,qintptr,qpackn,qvalue,
     :                 rel,skip)

      if (.not. skip) then
         call lsp_comp(qnoc,qj1,iqlsp,index,nterm,
     :               njv,maxj,minj,ncfg,iuc,pflsj,plnze)
      end if

      nze = 0
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      do jb = myid+1,ncfg,nprocs
        call SHELLSJB(JB)
        call BreitGG(NEW,NZERO,IFIRST,idg,irfst,skip,nze,plnze)
        if (mod(jb,100) .eq. 0) write(0,*) '   jb =',jb
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end do

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
291   NWFD = MAXORB
      if (qiajcmp.ne.0) call dalloc(qiajcmp,nwfd)
      if (qljcomp.ne.0) call dalloc(qljcomp,nwfd)
      if (qnjcomp.ne.0) call dalloc(qnjcomp,nwfd)
      if (qiajcld.ne.0) call dalloc(qiajcld,nwcd)
      if (qljclsd.ne.0) call dalloc(qljclsd,nwcd)
      if (qsignfa.ne.0) call dalloc(qsignfa,ncfg)
      if (qirfst.ne.0) call dalloc(qirfst,ncfg)
      if (qacmult.ne.0) call dalloc(qacmult,nwfd)
      if (qnelcsh.ne.0) call dalloc(qnelcsh,8*ncfg)
      if (qnocorb.ne.0) call dalloc(qnocorb,8*ncfg)
      if (qnoc.ne.0) call dalloc(qnoc,ncfg)
      if (qj1.ne.0) call dalloc(qj1,15*ncfg)
*
      return
      END
