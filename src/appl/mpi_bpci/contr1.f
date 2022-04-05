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
         if (idisk .eq. 1) rewind(iouhm)
         do 30 j = myid+1,ncfg,nprocs
            if (idisk .eq. 1) then
	       nij = 0
	    end if
	    if (mod(lsp(j),2) .eq. 0) then
	       nij = nij+1
	       else
	          read(iouhn) jb
	       end if
 	    else
	       if (iipos.le. ncfg) then
	          nij = nij + 1
	       end if
	    end if
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
             do j = myid+1,ncfg,nprocs
	        mycol = mycol + 1
	        read(iouhn) jb,m1,(hmx(nij+i),i=1,m1),(imx(nij+i),i=1,m1)
	        hii(j) = hmx(nij+1)
	        nij = nij + m1
	        jptr(mycol) = nij
             end do
	     do j = 1+myid+nprocs,ncfg,nprocs
                mycol = mycol + 1
                hmx(jptr(mycol-1)+1) = hmx(jptr(mycol-1)+1)-shift
             end do
	  end if
      end if
**********************end iLS = 1 *********************

