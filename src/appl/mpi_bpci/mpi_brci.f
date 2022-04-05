*     ==================================================================
*       A CONFIGURATION INTERACTION PROGRAM EITHER NON-RELATIVISTIC
*             OR IN THE BREIT-PAULI APPROXIMATION
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 1983
*
*       Modified August, 1992 to combine BREIT and the sparse
*       matrix Dvdson eigensolver.
*** Modified for MPI May 1996 , by Andreas Stathopoulos
*
*     ==================================================================
*
*       The PARAMETER values in this program define the following:
*               NOD   - Maximum number of points in the range of a
*                         - function
*               LIMD  - Maximum number of basis vectors
*               LCDIM - Number of array elements in a block
*                       of the direct access file
*
*       Variables for these dimensions are
*               NWF   - Number of functions (or electrons)
*               NO    - Number of points in the range of the function
*               NCFG  - Number of configuration state functions
*               NUME  - Number of eigenvalues/eigenvectors
*               NZE   - Number of non-zero matrix elements
*
*     ------------------------------------------------------------------
*      PROGRAM BRCI
      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60, NTERMD=31)
* 
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
        Character*2 idstring
        character*128 NAME(2)
****************************************************************
*
      character*24 strng
      character*128 ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      equivalence (ATOMC,NAME(1)),(ATOML,NAME(2))
      integer ibogus(16)
      equivalence (ATOMC,ibogus)

      CHARACTER QREL, QMASS
      character str*72
      LOGICAL REL, SKIP
      logical jhere,lhere,mhere,onlydvd
      logical leigen(NTERMD,NTERMD)
      integer nume(NTERMD)
      logical inputOK
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER NJEIG(ntermd)
      POINTER (IQLSP,LSP(1))
      pointer(pflsj,nze_flsj(1,1,1,1));
      double precision nze_flsj
* Start mpi!
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*     Set machine dependent unit numbers
*               IREAD - The unit number of standard input
*               IWRITE- The unit number of standard output
*               ISCW  - The unit number of standard error (screen)
*               iuc   - The unit number for the <name>.c file
*               iuw   - The unit number for the <name>.w file
*               ioul  - The unit number for the <name>.l file
*               iouj  - The unit number for the <name>.j file
*               iouhn  - The unit number for storing the LS matrix
*               iouhz - The unit number for storing the ZNV data
*               iouhs - The unit number for storing the S data
*               iouhm - The unit number for writing the scratch matrix
*
      write(idstring,'(I2.2)') myid
      if (myid.eq.0) then
        ISCW = 0
        IWRITE = 6
      else
        ISCW = 30
        IWRITE = 30
      endif
       !call alloc(iqlsp,30,4)
      IREAD = 5
      iuc    = 4
      iuw    = 2
      ioul   = 7
      iouj   = 8
      iounew = 13
      iouhn  = 9
      iouhz  = 10
      iouhs  = 11
      iouhm  = 12

 9    FORMAT(//20X,'==================================='/
     :         25X,' BRCI_MPI ... 2000   '/
     :         20X,'==================================='/)
      write(iscw,9)
      write(iscw,*) '                 ...brci_mpi running on ',
     :                  nprocs,' processors...'
      write(iscw,*)
 
      CALL INITA
      CALL INITR
      REL = .true.
* Node 0 executes the following IO * * * * * * * * * * * * * * * * * * *
*
      if (myid .eq. 0) then

1     WRITE(ISCW,'(/A,A)') ' Enter ATOM, relativistic (Y/N)',
     :               ' with mass correction (Y/N)'
      READ(IREAD,'(A)') STRNG
      I = INDEX(STRNG,',')
      IF (I .eq. 0)     THEN
         WRITE(ISCW,*) ' Separate with commas'
         GO TO 1
      ELSE IF (I .GT. 21) THEN
         WRITE(ISCW,*) ' ATOM name can have at most 20 characters'
         GO TO 1
      END IF
      QREL  = STRNG(I+1:I+1)
      QMASS = STRNG(I+3:I+3)
      ATOMC = STRNG(1:I-1)//'.c'
      ATOML = STRNG(1:I-1)//'.l'
      ATOMJ = STRNG(1:I-1)//'.j'
      ATOMW = STRNG(1:I-1)//'.w'
      ATOMNEW = STRNG(1:I-1)//'.new'
      IF (QREL.EQ.'N' .OR. QREL.EQ.'n') REL = .FALSE.
      MASS = 0
      IF (QMASS.EQ.'Y' .OR. QMASS.EQ.'y') THEN
         WRITE(ISCW,'(A)') ' Gradient or Slater form? (G/S):'
         READ(IREAD,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      END IF

*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	WRITE(ISCW,*) ' Use existing Matrix and ',
     :			'<atom>.l/j initial guess (Y/y)?'
	READ(IREAD,'(A)') STRNG
	if (STRNG(1:1).eq.'Y'.or.STRNG(1:1).eq.'y') then 
           onlydvd = .true.
	   inquire(file=ATOML,EXIST=lhere)
	   inquire(file=ATOMJ,EXIST=jhere)
           OPEN(UNIT=iounew,FILE=ATOMNEW,STATUS='UNKNOWN')
	   write(ISCW,*) ' The file <>.new will contain',
     :		' the new <>.l or <>.j info.'
	endif
*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*    Units iuc and iuw will be CLOSED by BREVAL
*
        OPEN(UNIT=iuw,FILE=ATOMW,STATUS='OLD', FORM='UNFORMATTED')
        OPEN(UNIT=ioul,FILE=ATOML,STATUS='UNKNOWN')
        OPEN(UNIT=iouj,FILE=ATOMJ,STATUS='UNKNOWN')
      ENDIF
*	..send the atomc string
      ilen = i+2
      call MPI_BCAST(ilen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(onlydvd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lhere,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jhere,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CPJ      call MPI_BCAST(atomc,ilen,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(atomc,128,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
c     to avoid MPI bug:
*      call MPI_BCAST(ibogus,ilen,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
	   
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* ALL NODES DO: everybody will read iuc (the clist)
      OPEN(UNIT=iuc,FILE=ATOMC,STATUS='OLD')

*	 * * * * * * Test if matrix exists * * * * * * * * * * *
	inquire(file='hnr.lst.'//idstring,EXIST=mhere)
	if (onlydvd) then 
           if (.not.mhere) call MPI_FINALIZE(ierr) 
	   if (jhere) then
	       inquire(file='hzeta.lst.'//idstring,EXIST=mhere)
	       if (.not.mhere) call MPI_FINALIZE(ierr)  
	       inquire(file='hspin.lst.'//idstring,EXIST=mhere)
	       if (.not.mhere) call MPI_FINALIZE(ierr) 
	   endif
	endif
*	 * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      OPEN(UNIT=iouhn,FILE='hnr.lst.'//idstring,STATUS='unknown',
     :        form='unformatted')
      OPEN(UNIT=iouhz,FILE='hzeta.lst.'//idstring,STATUS='unknown',
     :        form='unformatted')
      OPEN(UNIT=iouhs,FILE='hspin.lst.'//idstring,STATUS='unknown',
     :        form='unformatted')

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Only node 0 does I/O of J values and eigen values
      if (myid .eq. 0) then
         inputOK = .false.
         do while (.not.inputOK) 
            WRITE(ISCW,*) ' Enter Maximum and minimum ',
     :                           'values of 2*J'
            READ(IREAD,*) MAXJ,MINJ
            njv = (maxj-minj)/2 + 1
            if (njv .gt. 20) then
              write(iscw, *) ' Current dimensions ', 
     :                    'allow upto 20 J-values: Re-enter'
            else if (maxj < minj) then
              write(iscw, *) ' max must be greater than max: Re-enter'
            else
              inputOK = .true.
            end if
         end do
         ij = 1
         maxnum = 1
         nb = 1
         leigen = .false.
         write(iscw,*)
         write(iscw,*) 'Enter eigenvalues: ',
     :          'one line per term, eigenvalues separated by commas'
         do nj = maxj,minj,-2
           write(iscw,'(A,I3)') '2*J =',nj
           read(iread,'(A)') str
           nch = 1
           len = len_trim(str)
           do while (nch <= len)
             ipos = index(str(nch:len),',')
             if (ipos .eq. 0) ipos = len+2-nch
             read (str(nch:nch+ipos-2),*) keigv
             if (keigv .gt. ntermd) then
               write(0,*) 'Too high an eigenvalue requested:',
     :                   'Maximum for current dimensions is',ntermd
               stop
             end if
             leigen(keigv,nb) = .true.
             nch = nch + ipos
           end do
           nume(nb) = keigv
           njeig(nb) = nume(nb)
           if (nume(nb).gt. maxnum) maxnum = nume(nb)
           nb = nb + 1
         end do
      endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

********************************************************************
*	.. send to other nodes
      call MPI_BCAST(rel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(minj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mass,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxnum,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(njeig,ntermd,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(leigen,NTERMD*NTERMD,MPI_LOGICAL,0,
     :                 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nume,NTERMD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(njv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
      CALL BREVAL(rel,skip,ec,nze,onlydvd,iounew,njv,maxj,minj,
     :            pflsj)

      print*, 'All configurations processed: MYID = ', myid
      if (.not.onlydvd) then
        endfile(iouhn)
        endfile(iouhz)
        endfile(iouhs)
      else
	nze = 0
	do i=myid+1,NCFG,nprocs
	   read(iouhn) jb,m
	   nze = nze+m
	enddo
      endif

      rewind(iouhn)
      rewind(iouhz)
      rewind(iouhs)

*     .. allocate the memory
!      call alcmat(ncfg,nze,maxnum)
      if (skip) then
         iLS = 1
         close(unit=iouhz,status='delete')
         close(unit=iouhs,status='delete')
!         if (idisk .eq. 1) iouhm = iouhn
      else
         iLS = 0
!         if (idisk .eq. 1)
!     :         open(iouhm,file='mat.tmp.'//idstring,
!     :              status='unknown',form='unformatted')
      end if
!             call EXIT(5)

      ij = 1
      do j = maxj, minj,-2
	call lsjmat(j,nume(ij),skip,nclosd,nwf,ncfg,nze,ec,
     :    onlydvd,jhere,lhere,iounew,leigen(1,ij),pflsj,njv,maxj)
	ij = ij+1
      end do
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	if (myid.eq.0) then
        write(iounew,*) '**'
	endif
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if (idisk .eq. 1) close(iouhm,status = 'delete')
      CALL  MPI_FINALIZE(ierr)
      END
