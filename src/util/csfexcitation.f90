!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CSFEXCITATION
!
! This program generates excitation input to lsgen
!
! Per Jönsson & Jorgen Ekman, Malmo University, August 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program csfexcitation
implicit none
integer            :: i,k,j,nmax,lmax,mr,n(15),l(15),occ(15),nr,lr,occr,nl(0:15)
integer            :: jmin,jmax,nexc,occupied,ncore,nstart(0:6),lstart(0:6)
integer            :: nc(30),lc(30),ncoreloop,nclose,norbitalstring
integer            :: conflength,norb,orblength(15),orbstart(0:15),flag(15),nfound
integer            :: number1,number2,norbstrings,orbstring_l(11),orbstring_n(11)
integer            :: ncoreorbitals,lmaxcore,ntimes,ndouble

integer            :: pos,jl,jr

character(len=1)   :: lstring,ans
character(len=2)   :: sel(15)
character(len=135) :: config
character(len=135) :: configvect(100)
character(len=44)  :: orbitalstring
character(len=11)  :: decoding
character(len=3)   :: orbital(11)
character(len=2)   :: term


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

decoding = 'spdfghiklmn'

open(11, file='excitationdata.sh',status='unknown',form='formatted')
open(12, file='csfexcitation.log',status='unknown',form='formatted')
call system('chmod 755 excitationdata.sh') 

! Generate initial input to lsgen

write(11,'(a)') '#!/bin/bash'

nc = 0
lc = 0
ncoreorbitals = 0
lmaxcore = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) 
write(*,*) 'CSFEXCITATION'
write(*,*) 'This program creates excitation input to LSGEN'
write(*,*) 'Configurations should be entered in spectroscopic notation'
write(*,*) 'with occupation numbers and indications if orbitals are'
write(*,*) 'closed (c), inactive (i), active (*) or has a minimal'
write(*,*) 'occupation e.g. 1s(2,1)2s(2,*)'
write(*,*) 'Outputfiles: excitationdata.sh, csfexcitation.log'
write(*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ntimes = 1
110 continue
n = 0
l = 0
occ = 0

! Input the mutireference for later processing

write(*,*) 'Enter list of (maximum 100) configurations. End list with a blank line or an astersik (*).'
write(*,*) 

do j = 1,100

10 continue
   write(*,*) 'Give configuration', j
   read(*,'(a)') configvect(j)
   if(trim(configvect(j)).eq.''.or.trim(configvect(j)).eq.'*') then
      mr = j-1
      write(12,'(a)') '*'
      goto 99
   end if
   write(12,'(a)') trim(configvect(j))

!  Initial check, each orbital need to be closed, inactive, or minimal

   conflength = len(trim(configvect(j)))
   number1 = 0
   number2 = 0
   do k = 1,conflength
      if (configvect(j)(k:k).eq.'(') number1 = number1 + 1 
      if (configvect(j)(k:k).eq.',') number2 = number2 + 1 
   end do
   if (number1.ne.number2) then
      write(*,*)  'Each orbital must be closed (c), inactive (i), active (*)'
      write(*,*)  'or have a minimal occupation; redo!'
      goto 10
   end if
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input and process orbital set

99 continue
write(*,*) 'Give set of active orbitals in a comma delimited list ordered by l-symmetry, e.g. 5s,4p,3d'
read(*,'(a)') orbitalstring
write(12,'(a)') trim(orbitalstring)

! De-code orbital string

! Checks on orbital string

jl=1;
jr=0;
k=1;
do while ( index(orbitalstring(jl:len_trim(orbitalstring)),',').gt.0 )
   pos = index(orbitalstring(jl:len_trim(orbitalstring)),',')
   jr = jr + pos
   if(len_trim(orbitalstring(jl:jr-1)).eq.3) then
      orbital(k) = orbitalstring(jl:jr-1)
   else if(len_trim(orbitalstring(jl:jr-1)).eq.2) then
      orbital(k) = " " // orbitalstring(jl:jr-1)
   else
      write(*,*) 'Orbitals should be given in comma delimited list, redo!'
      goto 99      
   end if
   jl = jr + 1
   k = k +1
end do

jr = len_trim(orbitalstring)
if(len_trim(orbitalstring(jl:jr)).eq.3) then
   orbital(k) = orbitalstring(jl:jr)
else if(len_trim(orbitalstring(jl:jr)).eq.2) then
   orbital(k) = " " // orbitalstring(jl:jr)
else
   write(*,*) 'Orbitals should be given in comma delimited list, redo!'
   goto 99      
end if

norbstrings = k

! Determine highest n for orbital set and highest n for each symmetry

nmax = 0
do i = 1,norbstrings
!   write(*,'(a)') orbital(i)
select case(orbital(i)(3:3))
   case('s')
      lmax = 0
   case('p')
      lmax = 1
   case('d')
      lmax = 2
   case('f')
      lmax = 3
   case('g')
      lmax = 4
    case('h')
      lmax = 5
   case('i')
      lmax = 6
   case('k')
      lmax = 7
   case('l')
      lmax = 8
   case('m')
      lmax = 9
   case('n')
      lmax = 10
   case default
      write(*,*) 'Orbital quantum numbers should be in range s to n'
      stop
end select
   orbstring_l(i) = lmax

select case(orbital(i)(1:2))
   case(' 1')
      nmax = 1
   case(' 2')
      nmax = 2
   case(' 3')
      nmax = 3
   case(' 4')
      nmax = 4
   case(' 5')
      nmax = 5
   case(' 6')
      nmax = 6
   case(' 7')
      nmax = 7
   case(' 8')
      nmax = 8
   case(' 9')
      nmax = 9
   case('10')
      nmax = 10
   case('11')
      nmax = 11
   case('12')
      nmax = 12
   case('13')
      nmax = 13
   case('14')
      nmax = 14
   case('15')
      nmax = 15
   case default
      write(*,*) 'Principal quantum numbers should be <= 15'
      stop
end select
   orbstring_n(i) = nmax
end do

! Determine highest n and l for orbital set as well as spectrocopic notation

nmax = maxval(orbstring_n(1:norbstrings))
if (lmax.lt.lmaxcore) lmax = lmaxcore
lstring = decoding(lmax+1:lmax+1)

do i = 0,lmax
   nfound = 0
   do j = 1,norbstrings
      if (orbstring_l(j).eq.i) then
        nl(i) = orbstring_n(j)
        nfound = 1
      end if
   end do
   if (nfound.eq.0) nl(i) = 0
end do

if (lmax.ge.nmax) then
   write(*,*) 'Orbital quantum number should be less than n'
   stop
end if

if(ntimes.eq.1) then
   write(*,*) 'Resulting term? (1S, 3P, etc.)'
   read(*,*) term
   write(12,'(a)') term
end if
write(*,*) 'Number of excitations (if negative number e.g. -2, correlation '
write(*,*) 'orbitals will always be doubly occupied)                        '
read(*,*) nexc
ndouble = 0
if (nexc.lt.0) then
   ndouble = 1
end if
write(12,*) nexc, ' ! Number of excitations '

!write(*,*)
!write(*,*) 'TESTING ORBITAL SET'
!do i = 0,lmax
!   write(*,*) 'orbital set, l, nmax',i,nl(i)
!end do
!write(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOOPSTART
do j = 1,mr
   write(11,'(a)') 'lsgen << EOF'

! New list or add to list
   if(j.eq.1.and.ntimes.eq.1) then
      write(11,'(a)') '*'
   else
      write(11,'(a)') 'a'
   end if

! MCHF
   write(11,'(a)') '*'

! Default symmetry
   write(11,'(a)') '*'

! Use new reference set. Only if j > 1 or ntimes > 1
   if(j.gt.1.or.ntimes.gt.1) then
      write(11,'(a)') '*'
   end if

! Highest principal quantum number
   write(11,*) nmax

! Highest orbital quantum number
   write(11,'(a)') lstring

! All these orbitals active
   write(11,'(a)') 'n'

! Limitations of population
   write(11,'(a)') '*'

   config = configvect(j)

! Analyze config and extract n, l and occupation

   ! Determine length of orbit substrings. Depends on number of positions
   ! that occupation and type indicator occupies.
   conflength = len(trim(config))
   k = 0
   orbstart(0) = 0
   do i = 1, conflength
      if(config(i:i).eq.')') then
         k = k + 1
         orblength(k) = i - orbstart(k-1)
         orbstart(k) = i
      end if
   end do
   norb = k  ! Number of orbit substrings in configuration

   flag(:) = 0
   do k = 1,norb
      i = 1+orbstart(k-1)
      if(config(i+4:i+4).eq.',') flag(k) = 1 ! Flag to facilitate determination if  
                                             ! occupation/type indicator occupies one/two positions
      select case(config(i:i))
         case('1')
            nr = 1
         case('2')
            nr = 2
         case('3')
            nr = 3
         case('4')
            nr = 4
         case('5')
            nr = 5
         case('6')
            nr = 6
         case('7')
            nr = 7
         case('8')
            nr = 8
         case('9')
            nr = 9
      end select

      select case(config(i+1:i+1))
         case('s')
            lr = 0
         case('p')
            lr = 1
         case('d')
            lr = 2
         case('f')
            lr = 3
         case('g')
            lr = 4
         case('h')
            lr = 5
         case('i')
            lr = 6
         case('k')
            lr = 7
         case('l')
            lr = 8
      end select

      ! If occupation number "occupies" one position
      if(flag(k).eq.1) then
         select case(config(i+3:i+3))
           case('1')
              occr = 1
           case('2')
              occr = 2
           case('3')
              occr = 3
           case('4')
              occr = 4
           case('5')
              occr = 5
           case('6')
              occr = 6
           case('7')
              occr = 7
           case('8')
              occr = 8
           case('9')
              occr = 9
          end select
      else
         select case(config(i+3:i+4))
           case(' 1')
              occr = 1
           case(' 2')
              occr = 2
           case(' 3')
              occr = 3
           case(' 4')
              occr = 4
           case(' 5')
              occr = 5
           case(' 6')
              occr = 6
           case(' 7')
              occr = 7
           case(' 8')
              occr = 8
           case(' 9')
              occr = 9
           case('10')
              occr = 10
           case('11')
              occr = 11
           case('12')
              occr = 12
           case('13')
              occr = 13
           case('14')
              occr = 14
          end select
       end if
      if (nr.gt.nmax) then
         write(*,*) 'n in config greater than nmax'
         stop
      end if

      n(k) = nr
      l(k) = lr
      occ(k) = occr

      ! Determine start position and number of positions for type indicator
      if(orblength(k).eq.7) then
         sel(k) = config(i+5:i+5)
      elseif((orblength(k).eq.8).and.(flag(k).eq.0)) then
         sel(k) = config(i+6:i+6)
      elseif((orblength(k).eq.8).and.(flag(k).eq.1)) then
         sel(k) = config(i+5:i+6)
      else
         sel(k) = config(i+6:i+7)
      end if

!CPJ  Orbitalerna i configurationen med occupation and indicators
!      write(*,*) 'n,l,occ',n(k),l(k),occ(k), sel(k)

   end do

! Generate input to jjgen from this configuration

! Highest n for reference configuration
   write(11,*) maxval(n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   write(*,*) 'TESTING CONFIG',j

!   do k = 1, ncoreorbitals
!      write(*,*) 'coreorbitals',k,nc(k),lc(k)
!   end do

!   do k = 1, norb
!      write(*,*) 'orbitals with occupation ',k,n(k),l(k),occ(k),sel(k)
!   end do

!   do nr = 1,maxval(n)
!      do lr = 0,min(nr-1,lmax)
!         write(*,*) 'orbitals upp to nmax of configuration',nr,lr
!      end do
!   end do

!   do nr = maxval(n)+1,nmax
!      do lr = 0,min(nr-1,lmax)
!         write(*,*) 'orbitals upp to total nmax',nr,lr
!      end do
!   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do nr = 1,maxval(n)
      do lr = 0,min(nr-1,lmax)

! Check if this orbital was previously defined as closed. If so do not write anything for this

          nclose = 0
          do k = 1,ncoreorbitals
             if ((nr.eq.nc(k)).and.(lr.eq.lc(k))) then
                nclose = 1
             end if
          end do
          if (nclose.eq.1) cycle

! Find out if orbital is closed then add it to core orbitals 
! In for the next reference configuration there will be no question for this

          do k = 1, norb
             if ((nr.eq.n(k)).and.(lr.eq.l(k))) then
                if (trim(adjustl(sel(k))).eq.'c') then 
                   ncoreorbitals = ncoreorbitals + 1
                   nc(ncoreorbitals) = n(k)
                   lc(ncoreorbitals) = l(k)
                end if
             end if
          end do

! Find out occupation and options for this orbital

          occupied = 0
          do k = 1, norb
             if ((nr.eq.n(k)).and.(lr.eq.l(k))) then
                write(11,*) occ(k)
                write(11,'(a)') trim(adjustl(sel(k)))
                occupied = 1
             end if
          end do
          if (occupied.eq.0) then
             write(11,*) 0
!CPJ             write(11,'(a)') '*'
!Correction
             if (nr.le.nl(lr)) then
                if (ndouble.eq.1) then
                   write(11,'(a)') 'd'
                else
                   write(11,'(a)') '*'   
                end if
             else
                write(11,'(a)') 'i'
             end if
!Correction
          end if

      end do
   end do

! Loop over remaining orbitals in the orbital set

   if (nmax.gt.maxval(n)) then
      do nr = maxval(n)+1,nmax
         do lr = 0,min(nr-1,lmax)
            if (nr.le.nl(lr)) then
               if (ndouble.eq.1) then
                  write(11,'(a)') 'd'    
               else
                  write(11,'(a)') '*'
               end if
            else
               write(11,'(a)') 'i'
            end if
         end do
      end do
   end if

! Write term  and number of excitations
   if(j.eq.1.and.ntimes.eq.1) then
      write(11,'(a)') term
   end if
   write(11,*) iabs(nexc)

   if (j.ne.mr) then
      write(11,'(a)') '*'
   else
      write(*,'(a)') ' Generate more lists ? (y/n)'
      read(*,*) ans
      write(12,'(a)') ans
      if ((ans.eq.'y').or.(ans.eq.'Y')) then
         ntimes = ntimes + 1
         write(11,'(a)') '*'
         write(11,'(a)') 'EOF'
         write(11,*)
         write(11,'(a)') 'cp clist.out clist.inp'
         write(11,*)
         write(11,'(a)') 'lsgen << EOF'
         write(11,'(a)') 'r'
         write(11,'(a)') 'EOF'
         write(11,*)
         write(11,'(a)') 'cp clist.out clist.inp'
         write(11,*)
         goto 110
      else
         write(11,'(a)') '*'
      end if
   end if
   write(11,'(a)') 'EOF'

   if(mr.gt.1.and.j.lt.mr) then
      write(11,*)
      write(11,'(a)') 'cp clist.out clist.inp'
      write(11,*)
      write(11,'(a)') 'lsgen << EOF'
      write(11,'(a)') 'r'
      write(11,'(a)') 'EOF'
      write(11,*)
      write(11,'(a)') 'cp clist.out clist.inp'
      write(11,*)
   end if
end do

close(11)


end program csfexcitation


