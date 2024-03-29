      subroutine mpi_hmx_lsj(ncfg,jj,j,nze,ind_a,nij,istart,shift,
     :    idisk,mycol,pflsj,maxj)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NTERMD=31,NWD=60)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************
      LOGICAL  SKIP
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,termsh(ntermd),nterm
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /DSSLSJ/ shift
*
      INTEGER cptr
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),(iqkval,kval(1)),(iqcptr,cptr(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
      INTEGER   q(8), ipos(3),iselec(NTERMD)
      LOGICAL iupper,hiend
      DATA iupper/.false./
      DATA crite,critc,critr,trhold,ortho/
     :       0.d-17,1.d-12,1.d-20,1.d0,1.d0/
      pointer(pflsj,flsj(nterm,nterm,2,maxj))
      double precision flsj;
      logical lhmx_memory,lhmx_disk

*
         jb_term = 0;
         do while (index(jb_term+1) <= j)
            jb_term = jb_term + 1;
         if (jb_term == nterm) exit;
         end do
         mycol = mycol + 1
         if (idisk .eq. 1) then
	    nij = 0
            call dfill(nze,0.d0,hmx,1)
	 end if
	 if (mod(lsp(jb_term),2) .eq. 0) then
*            .. csf does not contribute
	    nij = nij+1
	    imx(nij) = j
	    jptr(mycol) = nij
	    hmx(nij) = 1.d+5
	    if (j .eq. 1) then
	       read(iouhn) jb,m1,(h(i,1),i=1,m1),
     :                           (jan(i,1),i=1,m1)
	    else
	       read(iouhn) jb
	    end if
	    read(iouhz) jb
	    read(iouhs) jb
 	 else
            ipos(1:3) = 0
	    ntermj = mod(lsp(jb_term),64)/2
	    read(iouhn) jb,m1,(h(i,1),i=1,m1),
     :                        (jan(i,1),i=1,m1)
	    read(iouhz) jb,m2,(h(i,2),i=1,m2),
     :                        (jan(i,2),i=1,m2)
            read(iouhs) jb,m3,(h(i,3),i=1,m3),
     :                        (jan(i,3),i=1,m3)
            call advance(1,m1,ipos,ncfg)
	    call advance(2,m2,ipos,ncfg)
	    call advance(3,m3,ipos,ncfg)
*           .. while there is more data
*             find which file has lowest row
            iipos = min(ipos(1),ipos(2),ipos(3))
            do while (iipos <= ncfg)
*             .. there is more data
	       nij = nij + 1
	       irow = min(nrow(1),nrow(2),nrow(3))
               ja_term = 0;
               do while (index(ja_term+1) <= irow)
                  ja_term = ja_term + 1;
                  if (ja_term == nterm) exit;
               end do
	       if (nrow(1) .eq. irow) then
	          hmx(nij) = hmx(nij) + h(ipos(1),1)
	          imx(nij) = jan(ipos(1),1)
	          call advance(1,m1,ipos,ncfg)
	       end if
	       if (nrow(2) .eq. irow) then
	          ntermi = mod(lsp(ja_term),64)/2
	          hmx(nij)= hmx(nij)+h(ipos(2),2)*
     :                          flsj(ntermi,ntermj,1,jj)
	          imx(nij) = jan(ipos(2),2)
	          call advance(2,m2,ipos,ncfg)
	       end if
	       if (nrow(3) .eq. irow) then
	          ntermi = mod(lsp(ja_term),64)/2
	          hmx(nij)= hmx(nij)+h(ipos(3),3)*
     :                          flsj(ntermi,ntermj,2,jj)
                  imx(nij) = jan(ipos(3),3)
	          call advance(3,m3,ipos,ncfg)
	       end if
               iipos = min(ipos(1),ipos(2),ipos(3));
            end do
	 end if
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         if (shift .eq. 0.d0) then
	    if (myid .eq. 0) then
               if (mod(lsp(jb_term),2) == 0) then
*             .. csf does not contribute
*             .. shift is determined from first element of hnr.lst
                  shift   = h(1,1)
	       else
                  shift = hmx(1)
	       end if
            endif
            call MPI_BCAST(shift,1,MPI_DOUBLE_PRECISION,0,
     $			MPI_COMM_WORLD,ierr)
	  endif

	  hmx(istart+1) = hmx(istart+1)-shift
	  hii(j) = hmx(istart+1)
	  if (idisk .eq. 1) then
	     write(iouhm) j,nij,(hmx(i),i=1,nij),
     :                          (imx(i),i=1,nij)
	  else
	     m1 = nij-istart
	     jptr(mycol) = nij
	     istart = nij
	  end if
     
      end subroutine hmx_lsj0
