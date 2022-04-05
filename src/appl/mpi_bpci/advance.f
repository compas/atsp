*
*     ------------------------------------------------------------------
*             A D V A N C E
*     ------------------------------------------------------------------
*
*     Advance the position of array j to the next relevant position
* 
      SUBROUTINE advance(j,len,ipos,ncfg)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NTERMD=31)
*
      DIMENSION ipos(3)
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,termsh(ntermd),nterm
      POINTER (IQLSP,LSP(1))
      pointer(pind,ind(nt)),(pnze_term,nze_term(njv)),
     :       (pnze_col,nze_col(njv)),(pnze_t,nze_t(njv)),
     :       (pnze_c,nze_c(njv))
      common/im_nze/pflsj,pnze_term,pnze_col,pind,nt,mxj,
     :              mnj,njv

*
*     ------------------------------------------------------------------
*
      i = ipos(j)
   10 i = i+1
      if ( i .gt. len) then
	ipos(j) = ncfg+1
	nrow(j) = ncfg+1
	return
      else
	irow = jan(i,j)
        ja_row = 1;
        do while (ind(ja_row) <= irow)
          ja_row = ja_row + 1;
          if (ja_row == nt) exit;
        end do
	if (mod(lsp(ja_row),2) .eq. 0) then
	  go to 10
	else
	  ipos(j) = i
	end if
      end if
      nrow(j) = jan(i,j)
!      print *, ' File, pos, row:',j,i,nrow(j),len
      end
