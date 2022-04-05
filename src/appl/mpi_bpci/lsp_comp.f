*-----------------------------------------------------------------------
*     terms 
*-----------------------------------------------------------------------

      subroutine lsp_comp(qnoc,qj1,iqlsp,index,nterm,
     :              njv,maxj,minj,ncfg,iuc,pflsj,plnze)

      parameter (MTERM=31, MEIG=31)
      character string*64
      character term_bl(mterm)*3, term_tmp*3,term_comp*3;
      integer iuc,noccsh,j1qnrd;
      integer index(MTERM),q(8),ls,m,lli,ksi,lt;
      integer llj,jb,ksj,nterm,nterms; 
!      integer llj,iuc,jb,ksj,nterm,nterms; 
      CHARACTER elc(8)*3, couple(15)*3
      character (len = 64)	:: config_tmp
      character (len = 66)      :: config(MEIG)
      integer		 	:: MAX_read, lstring,iounew 
      pointer(pflsj,nze_flsj(1));
      pointer(qnoc,noccsh(1)), (QJ1,J1QNRD(15,1))
      double precision nze_flsj;
      logical lnew_term, lnze 
      POINTER (IQLSP,LSP(1)), (plnze,lnze(4,njv))
      integer lsp
      integer mmj

      !call alloc(iqlsp,MTERM,4)
      nterm = 0;
      lnew_term = .true.
*     ... read 2 blank lines from file 'clist'
      rewind (iuc);
      Read (iuc,'(A64)') string
      Read (iuc,'(A64)') string

      i = 1; index(1) = 1;
      term_tmp = '000';
      term_bl = '';
      do ir = 1, ncfg
        READ(iuc,'(8(1X,A3,1X,I2,1X))') (ELC(j),Q(j),j=1,8)
        READ(iuc,'(15(1X,A3))') (COUPLE(J),J=1,15)
        nocc = 0
        do while(ELC(nocc+1).ne. '   ')
           NOCC = NOCC + 1
        end do
        call pack(nocc,elc,q,couple,config_tmp)
        lstring = Len_Trim(config_tmp)

        if (config_tmp(lstring:lstring) == '0' ) then
           term_comp = config_tmp((lstring-2):(lstring-1));
        else 
           term_comp = config_tmp((lstring-1):lstring);
        end if
        if (ir == 1) term_bl(1) = term_comp;

        icomp = 1;
        do while(icomp <= i) 
          if (term_bl(icomp).eq.term_comp) exit
          icomp = icomp + 1;
        end do

        if(icomp > i) then
          term_bl(icomp) = term_comp;
          nterms = icomp;
          i = icomp
        end if;
        index(i) = ir;
      end do
       
      call alloc(plnze,nterms*njv,8)
      lnze(1:nterms,1:njv) = .false.

      call alloc(iqlsp,nterms,4)
      lsp(1:nterms) = 0;
      nterm = 0;
      do i = 1,nterms
        m = noccsh(index(i))
        ls = j1qnrd(2*m-1,index(i))/64
        lli = mod(ls,64) -1
        ksi = (ls/64) - 1
        lsp(i) = (64*lli+ksi)*64 + 2*i
        jjc = 1
        do jj = maxj,minj,-2
           if((jj>=abs(lli-ksi)).and.(jj<=(lli+ksi))) 
     :             lnze(i,jjc) = .true.        
           jjc = jjc + 1
        end do 
        nterm = nterm + 1;
      end do

      call alloc(pflsj,nterm*nterm*2*njv,8);
      nze_flsj(1:(nterm*nterm*2*njv)) = 0.0;
      call flsj_comp(ncfg,nterm,index,njv,maxj,minj,
     :                pflsj,lsp)
      return
      end

*   compute the flsj  aray:
      subroutine flsj_comp(ncfg,nterm,index,njv,maxj,minj,
     :               pflsj,lsp)
      double precision w1,w2,nze_flsj,flsj;
      integer nterm
      pointer(pflsj,nze_flsj(nterm,nterm,2,njv));
      integer ii,i,ls,ksi,li,jjj,ksj,llj,jj,lsp(nterm),
     :        index(nterm);
      integer jc,nt,njv,maxj,minj,nze_term,nze_col,tmp
      pointer(pflsj_save,flsj(1))
      pointer(pind,ind(nterm)),(pnze_term,nze_term(njv)),
     :       (pnze_col,nze_col(njv))
      common/im_nze/pflsj_save,pnze_term,pnze_col,pind,nt,mxj,
     :              mnj,nj
      integer ind_jj
  
      pflsj_save = pflsj;
      nt = nterm;mxj=maxj;mnj=minj;nj=njv;
      call alloc(pind,nterm,4)
      call alloc(pnze_term,nterm,4);
      call alloc(pnze_col,nterm,4);
      ind(1:nterm) = index(1:nterm);
      nze_term(1:njv) = ncfg;
      nze_col(1:njv)= 0;
      ksi = 0; lli = 0; ksj = 0; llj = 0;
      ind_jj = 0
      do jj = maxj,minj,-2
         ind_jj = ind_jj + 1
         DO II = 1,NTERM
            call lsp_cv(lsp(ii),lli,ksi);
            DO JJJ = 1,NTERM
               call lsp_cv(lsp(jjj),llj,ksj);
               PHASE = (-1)**((LLI+KSJ-JJ+LLI+KSI-JJ+LLJ+KSJ-JJ)/2) 
               CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,2,W1) 
               CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,4,W2);
               nze_FLSJ(II,JJJ,1,ind_jj) = PHASE*W1;
               nze_FLSJ(II,JJJ,2,ind_jj) = PHASE*W2
            end do
         end do
      end do
      return
      end

      subroutine lsp_cv(ti,lli,ksi)
      integer ti,lli,kli
        ls = ti/64;
        ksi = mod(ls,64);
        lli = (ls/64)
      return
      end 


