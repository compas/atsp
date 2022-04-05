      SUBROUTINE LSJMAT(JJ, nume,skip,nclosd,nwf,ncfg,nze,ec,
     :    onlydvd,jhere,lhere,iounew,leigen,pflsj,njv,maxj)

      if (iLS .eq. 0) then
         do i = 1,ncfg
	    ls = lsp(i)/64
	    ksi = mod(ls,64)
	    lli = ls/64
	    if (jj .ge. ABS(LLi-KSi) .and. jj .le. lli+ksi) then
               lsp(i) = 2*(lsp(i)/2) + 1
	    else
               lsp(i) = 2*(lsp(i)/2)
            end if
         end do
         nij = 0
         istart = 0
         shift = 0.d0
         rewind(iouhn)
         rewind(iouhz)
         rewind(iouhs)
         if (idisk .eq. 1) rewind(iouhm)
         mycol = 0
         do 30 j = myid+1,ncfg,nprocs
            mycol = mycol + 1
            if (idisk .eq. 1) then
	       nij = 0
               call dfill(nze,0.d0,hmx,1)
	    end if
	    if (mod(lsp(j),2) .eq. 0) then
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
               ipos(1) = 0
               ipos(2) = 0
               ipos(3) = 0
	       ntermj = mod(lsp(j),64)/2
	       read(iouhn) jb,m1,(h(i,1),i=1,m1),
     :                          (jan(i,1),i=1,m1)
	       read(iouhz) jb,m2,(h(i,2),i=1,m2),
     :                        (jan(i,2),i=1,m2)
               read(iouhs) jb,m3,(h(i,3),i=1,m3),
     :                        (jan(i,3),i=1,m3)
               call advance(1,m1,ipos,ncfg)
	       call advance(2,m2,ipos,ncfg)
	       call advance(3,m3,ipos,ncfg)
  100          iipos = min(ipos(1),ipos(2),ipos(3))
	       if (iipos.le. ncfg) then
	          nij = nij + 1
	          irow = min(nrow(1),nrow(2),nrow(3))
	          if (nij .gt. nze) then
	             newlen = nze + nze/4
	             call realloc(iqhmx,nze,newlen,8)
                     call realloc(iqimx,nze,newlen,4)
                     call dfill(newlen-nze,0.d0,hmx(nze+1),1)
	             nze = newlen
                  end if
	          if (nrow(1) .eq. irow) then
	             hmx(nij) = hmx(nij) + h(ipos(1),1)
	             imx(nij) = jan(ipos(1),1)
	             call advance(1,m1,ipos,ncfg)
	          end if
	          if (nrow(2) .eq. irow) then
	             ntermi = mod(lsp(irow),64)/2
	             hmx(nij)= hmx(nij)+h(ipos(2),2)*
     :                             flsj(ntermi,ntermj,1,jj)
	             imx(nij) = jan(ipos(2),2)
	             call advance(2,m2,ipos,ncfg)
	          end if
	          if (nrow(3) .eq. irow) then
	             ntermi = mod(lsp(irow),64)/2
 	            hmx(nij)= hmx(nij)+h(ipos(3),3)*
     :                             flsj(ntermi,ntermj,2,jj)
                     imx(nij) = jan(ipos(3),3)
	             call advance(3,m3,ipos,ncfg)
	          end if
	          go to 100
	       end if
	    end if
            if (shift .eq. 0.d0) then
	       if (myid .eq. 0) then
	          if (mod(lsp(j),2) .eq. 0) then
                     shift   = h(1,1)
	          else
                     shift = hmx(1)
	          end if
               endif
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
 30       continue
	  call mpi_allr_dp(hii,ncfg,tm)
***************iLS = 0***************************************
      else
	 if (idisk .eq. 1 ) then
             rewind(iouhn)
	     mycol = 0 
             do j = myid+1,ncfg,nprocs
	        mycol = mycol + 1
                ipos(1) = 0
     	        read(iouhn) jb,m1,hii(j)
             end do
	     call mpi_allr_dp(hii,ncfg,tm)
             shift = hii(1)
             do j = 1,ncfg
     	        hii(j) = hii(j) - shift
             end do
	  else
             nij = 0
	     mycol = 0
             rewind(iouhn)
             do j = myid+1,ncfg,nprocs
	        mycol = mycol + 1
	        read(iouhn) jb,m1,(hmx(nij+i),i=1,m1),(imx(nij+i),i=1,m1)
	        hii(j) = hmx(nij+1)
	        nij = nij + m1
	        jptr(mycol) = nij
             end do
             call mpi_allr_dp(hii,ncfg,tm)
	     shift = hii(1)
	     do j = 1,ncfg
                hii(j) = hii(j) - shift
             end do
             hmx(1) = hii(myid+1)
             mycol = 1
	     do j = 1+myid+nprocs,ncfg,nprocs
                mycol = mycol + 1
                hmx(jptr(mycol-1)+1) = hmx(jptr(mycol-1)+1)-shift
             end do
	  end if
      end if
**********************end iLS = 1 *********************

      nloc = mycol
      if (iLS .eq. 1) then
	iouv = ioul
	jhere = lhere
      else
	iouv = iouj
	lhere = jhere
      end if

      lim = min(20,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      iiwsz = 6*lim + nume
      n = ncfg
      ilow = -1
      ihigh = -1
      ie = 1

      do ic1 = 1,ntermd
         if (leigen(ic1)) then
           iselec(ie) = ic1
           ie = ie + 1
         end if
      end do

      iselec(ie) = -1
      maxiter = max(60*nume,150)
      mblock = ie -1
      niv = 0

      if (onlydvd .and. lhere) then
         if (lhere) then
            if (myid.eq.0) then
              rewind(iouv)
              read(iouv,*)
              read(iouv,'(//8X,I4,2X,8X,I4)') jj1,number
              niv = number
	      do K=0,number-1
	         read(iouv,*)
	         read(iouv,*)
	         read(iouv,'(7F11.8)') (EIGVEC(J+K*N),J=1,N)
	      enddo
	      kcur = 1
	      do k=1,number
                jcur = kcur + N
                call dgemv('T',N,number-k,1.d0,eigvec(jcur),N,
     :                  eigvec(kcur),1, 0.d0,tm,1)
                do j=k+1,number
                   call daxpy(N,-tm(j-k),eigvec(kcur),1,eigvec(jcur),1)
                   jcur = jcur + N
                enddo
                kcur = kcur + N
              enddo
*	      ..Normalize		- - - - - - - - - - - - 
	      kcur = 1
              do k=1,number
                tm(k) = 1./dnrm2(N,eigvec(kcur),1)
                call dscal(n,tm(k),eigvec(kcur),1)
                kcur = kcur + N
 	      enddo
	    endif
	    call MPI_BCAST(niv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(eigvec,number*N,MPI_DOUBLE_PRECISION,0,
     $			MPI_COMM_WORLD,ierr)
         endif
      endif

      CALL dvdson(hmx,imx,jptr,iupper,nze,tm,tp,
     :     Ncfg,lim,hii,ilow,ihigh,iselec,niv,mblock,
     :     crite,critc,critr,ortho,maxiter,
     :     eigvec,iworksz,iwork,iiwsz,hiend,nloops,nmv,ierr)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      END 
