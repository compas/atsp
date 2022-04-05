 
*     ------------------------------------------------------------------
*       R E A D W
*     ------------------------------------------------------------------
*
      SUBROUTINE READW(rel,ec,nnel,skip,onlydvd,iounew)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NWD=60,NOD=220)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************
	logical onlydvd
*     .. breit commons
      COMMON /CLOSED/B1ELC(4),NCLSD,IBK
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCF
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
*     .. radial commons
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      CHARACTER*3 el(NWD),el1
      CHARACTER*6 atom, term
      INTEGER iatom,iterm
c to avoid MPI bug
      EQUIVALENCE (atom,iatom),(term,iterm)
      LOGICAL found(NWD),check,rel,skip
      DOUBLE PRECISION PT(NOD)

      nwf = nclosd + maxorb
      call initr
*     copy from breit common to radial common
      do i = 1,nclosd
	 l(i) = ljclsd(i)
	 write(el(i),'(A3)') iajcld(i)
	 found(i) = .false.
      end do 
      do i = nclosd+1,nwf
	 l(i) = ljcomp(i-nclosd)
	 n(i) = njcomp(i-nclosd)
	 write(el(i),'(A3)') iajcmp(i-nclosd)
	 found(i) = .false.
      end do 
      ncfg = ncf
*
*  *****  READ THE RADIAL FUNCTIONS
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if (myid.eq.0) then
         IW = 1
  14     READ(iuw,END=13) ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
         call eptr(el,el1,i,*999)
         if (.not. found(i)) then
           do j= 1,m
              p(j,i) = pt(J)
           end do 
	   az(i) = azd
	   max(i) = m
           DO J = M+1,NO
              P(J,I) = D0
           end do
	   found(i) = .true.
           IW = IW+1
         end if
         Z = ZZ
         IF (IW.LE.NWF) GO TO 14
   13    close(iuw)
*     check that all functions were found
         check = .true.
         do i = 1,nwf
	     check = check .and. found(i)
         end do
         if (.not. check) then
	    write(iscw,*) ' Not all radial functions were found'
	    write(iscw,'(1X,A3,2X,L1)') (el(i),found(i),i=1,nwf)
	    call MPI_FINALIZE(ierr) 
         end if
      endif

	call MPI_BCAST(P,nod*nwf,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,ierr)
	call MPI_BCAST(AZ,nwf,MPI_DOUBLE_PRECISION,0,
     $		       MPI_COMM_WORLD,ierr)
	call MPI_BCAST(L,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(MAX,nwf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*	   ..send YK,YR,X from the common
	call MPI_BCAST(YK,3*nod,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,ierr)
*	   .. Z and TOL
	call MPI_BCAST(Z,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(iatom,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(iterm,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     .. initialize R and Rydberg constant
      call initm(1)
*     .. compute energy of core
      call ecore(el,ec,rel)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      if (myid.eq.0) then
         if (.not.onlydvd) close(unit=iounew, status='delete')
         if (skip) then 
	 if (.not.onlydvd) iounew = ioul
	 close(unit=iouj, status='delete')
      else
	 if (.not.onlydvd) iounew = iouj
	 close(unit=ioul, status='delete')
      end if

      write(iounew,'(2X,A6,A,F5.1,A,I3,A,I6)' ) ATOM,'  Z = ',Z,
     :        '  NEL = ', NNEL, '   NCFG = ',NCFG
	endif
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      return

  999 write(iscw,*) ' Electron',el1,' not found in electron list'
      write(iscw,*) ' List is',  el
      go to 14
      end
