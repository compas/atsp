*****:******************************************************************
*
* --- This SUBROUTINE allocates memory for the Interaction 
*     matrix
*
*****:******************************************************************
      Subroutine alcmat(skip,ncfg,nume,n_maxcol,n_hmx,lhmx_memory,
     :                  lhmx_disk,jj)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
* MINOR mpi MODIFICATIONS ************************************
*        
	INCLUDE 'mpif.h'
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),
     :        (iqkval,dummy1(1)),(iqcptr,dummy2(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      logical skip, lhmx_disk, lhmx_memory
      character*2 idstring
      integer n_maxcol, n_hmx
*
      write(idstring,'(I2.2)') myid
      lhmx_disk = .true.
      lhmx_memory = .true.
      idisk = 0;
      write (iscw,*)
      if (myid == 0) write (iscw,*) ':: Allocating memory for ',
     :                              'Block 2J = ',jj 
      lim = min(20,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      call alloc(iqeigvec,iworksz,8)
      iiwsz =  6*lim + nume
      call alloc(qtm,ncfg,8)
      call alloc(qtp,ncfg,8)
      call alloc(qdiag,ncfg,8)
      call alloc(qiwork,iiwsz,4)

      if (qtm.eq.0 .or. qtp.eq.0 .or. qdiag.eq.0 .or. qiwork.eq.0) then
         if (myid == 0) then
	   write(ISCW,*) 'Insufficient memory for diagonalization ',
     :                   ' failed at ##1(qtm,qtp,qdiag,qiwork)'
         end if
	 call MPI_FINALIZE(ierr)
      endif

      call alloc(iqjptr,ncfg,4)

      call alloc(iqhmx,n_hmx,8)
      call alloc(iqimx,n_hmx,4)
      if (iqhmx.eq.0.or.iqimx.eq.0) then
        lhmx_memory = .false.;
        if (myid == 0) then
          write(iscw,*) 'BLOCK 2J = ' , jj, 'Insufficient memory ',
     :          ' diagonalization ON DISK (hmx has ',n_hmx,') elements.'
        end if
      endif

      call MPI_ALLREDUCE(lhmx_memory,lhmx_memory,1,MPI_LOGICAL,
     :                     MPI_LAND, MPI_COMM_WORLD,ierror)

      if (.not.lhmx_memory) then
        if (iqhmx .ne. 0) call dalloc(iqhmx,n_hmx)
        if (iqimx .ne. 0) call dalloc(iqimx,n_hmx)
      end if

      if (.not.lhmx_memory) then
         idisk = 1;
         call alloc(iqhmx,n_maxcol,8)
         call alloc(iqimx,n_maxcol,4)
         if (iqhmx.eq.0.or.iqimx.eq.0) then
           lhmx_disk = .false.;
           if (iqhmx.ne.0) call dalloc(iqhmx,n_hmx)
           if (iqimx.ne.0) call dalloc(iqimx,n_hmx)
           if (myid == 0) then
             write(iscw,*) 'Insufficient memory for on-disk',
     :                     ' diagonalization: Exitting..'
           end if
           call MPI_FINALIZE(ierr)
         endif
      end if

      if (myid == 0) then
         if(lhmx_memory) then
             write(iscw,*) ':: IN MEMORY: Block 2J = ',jj, 'with ',
     :                      n_hmx, ' matrix', ' elements'
         end if
      end if

      if (skip) then
          if (idisk .eq. 1) iouhm = iouhn
      else
          if (idisk .eq. 1)
     :         open(iouhm,file='mat.tmp.'//idstring,
     :              status='unknown',form='unformatted')
      end if

      end subroutine alcmat



