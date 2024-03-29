*-----------------------------------------------------------------------
*      E I G _ O U T
*-----------------------------------------------------------------------

*     This routine writes first the eigenvectors, and then
*   the eigenvalues to Term.l 

      subroutine brci_eig(ncfg,iblock,nume,leigen,wt,Ssms,g_J,g_JLS,
     :                     eigval,EC,shift,iounew,iuc)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220,NOFFD=800)
      parameter (MTERM=31, MEIG=31)
      character string*64
      CHARACTER EL*3,ATOM*6,TERMs(MTERM)*6
      COMMON/LABEL/ EL(NWD),ATOM,TERMs
      logical leigen(nume)
      character term_bl(mterm)*3
      integer q(8),iounew,iuc
      CHARACTER elc(8)*3, couple(15)*3
      character (len = 64)	:: config_tmp
      character (len = 66), allocatable, dimension(:) :: config0
      character (len = 66)      :: config(MEIG)
      integer, dimension(2)	:: max_comp(MEIG)
      integer, dimension(1)	:: mx
      integer		 	:: MAX_read, lstring!,iounew 
      double precision wt(nume*ncfg),Ssms(nume),g_J(nume)
      double precision g_JLS(nume),EC,shift,eigval(nume)

      neig = 0
      max_comp = 0
      com_order = 0
      do i = 1,nume
         if (leigen(i)) neig = neig+1
      end do

*     ... read 2 blank lines from file 'clist'
      rewind (iuc);
      Read (iuc,'(A64)') string
      Read (iuc,'(A64)') string

*<<< find the configurations with max comp of the eigen vector
*     ... find the locations of extrema for all eigenvalues
      ics = 1
      ice = ncfg
      do ic1 = 1, nume
        eigval(ic1) = eigval(ic1) + EC + shift;
	mx = MaxLoc(abs(wt(ics:ice)))
        max_comp(ic1) = mx(1)
	ics = ics + ncfg
        ice = ice + ncfg 
      end do
*     ... find the maximum of max_comp to read max_read lines  
      max_read = MaxVal(max_comp)
*
*     allocate enough space to read the cfg
      allocate(config0(max_read));

*     ... read all lines.le.max_read of file 'clist'
      do i = 1, max_read 
        READ(iuc,'(8(1X,A3,1X,I2,1X))') (ELC(j),Q(j),j=1,8)
        READ(iuc,'(15(1X,A3))') (COUPLE(J),J=1,15)
        
*     ... find NOCC
        nocc = 0
        do while(ELC(nocc+1).ne. '   ')
           NOCC = NOCC + 1
        end do

*     ... pack config 
        call pack(nocc,elc,q,couple,config_tmp)
*        lstring = Len_Trim(config_tmp)
        config0(i) = trim(config_tmp)
      end do

*     ... order each config to correspond to max comp of the eig vector
      do i = 1, nume
        config(i) = config0(max_comp(i))
      end do
      
      deallocate(config0);
      jj = 0
      nnel = 0
*     ... write $TERM.l file for each block
      do i = 1, nume;
         if (leigen(i)) then;
           config_tmp = config(i) 
           lstring = Len_Trim(config_tmp)
           write (iounew,'(3(A8,f15.10))') 
     :         'Ssms=', Ssms(i), 'g_J=', g_J(i), 'g_JLS=', g_JLS(i)
           WRITE (iounew,'(i6,f16.9,2x,A)') 
     :                 max_comp(i), eigval(i), config_tmp(1:lstring) 
           WRITE (iounew,'(7F11.8)') (wt((i-1)*ncfg+j), j=1,ncfg)
*          WRITE(iscw,'(/1X,F19.12,2X,A/(1X,7F11.8))' ) EIGVAL
         end if
      end do

      return
      end
