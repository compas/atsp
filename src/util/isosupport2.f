      program isosupport
* 
* written by Per Jonsson and Michel Godefroid
* version July 31, 1998
*
      implicit double precision(a-h,o-z)

      parameter (emass = 548.579903d-6, ry = 10973731.534d0)
      parameter (a0 = 5.29177249d-11, pi = 3.141592654d0)
      doubleprecision isomass(10),el(10),eu(10),diff(10),eln(10),eun(10)
     :,rme(10),muovm(10),rydm(10)
      doubleprecision wlength(10),dwlength(10),rms(10)
      doubleprecision elns(10),euns(10)
      character*4 isolabel(10)
      character*72 transition

      open(unit=13,file='isodata',status='unknown')

      write(*,*) 
      write(*,*) '==================================================='
      write(*,*) ' Program for calculating transition isotope shift.'
      write(*,*) '==================================================='
      write(*,*) 
      write(*,*) '                                           _______'
      write(*,*) ' Observe that the root mean square radius \/ <r^2> '
      write(*,*) ' of the nucleus should be given, as separated from '
      write(*,*) ' the equivalent uniform nuclear radius Req where   '
      write(*,*) ' <r^2> = (3/5)*Req^2 '
      write(*,*) 
      write(*,*) ' Note that the elctron density |wfn(o)|^2 should '
      write(*,*) ' be given, as separated from the modified density '
      write(*,*) ' D = 4*pi*|wfn(0)|^2'
      write(*,*) '==================================================='

      write(13,*)
      write(13,*) '=============='
      write(13,*) ' Initial data'
      write(13,*) '=============='

      write(*,*) 
      write(*,*) 
      write(*,*) 
      write(*,*) ' Specify the transition'
      read(*,'(a72)') transition

      write(13,*)
      write(13,*) ' Transition'
      write(13,'(a72)') transition

      write(*,*) ' Nuclear charge'
      read(*,*) z

      write(13,*)
      write(13,*) ' Nuclear charge',z

      write(13,*)
      write(*,*) ' Energy (a.u.) infinite nucl. mass for lower state'
      read(*,*) elow
      write(13,*) ' Energy (a.u.) infinite nucl. mass for lower state',
     1 elow,'H'

991   write(*,*) 
     1' Energy (a.u.) infinite nucl. mass for upper state (1) or '
      write(*,*) 
     2' energy diff. (cm^-1) between upper and lower state (2)'
      read(*,*) n
      if (n.eq.1) then
        write(*,*) ' Energy (a.u.) infinite nucl. mass for upper state'
        read(*,*) eup
      elseif (n.eq.2) then
        write(*,*) ' Energy diff. (cm^-1) between upper and lower state'
        read(*,*) ediff
        eup = elow + 50.d0*ediff/ry
      else
         goto 991
      endif
      if (n.eq.1) then
         write(13,*) ' Energy for upper state from calculation'
      else
         write(13,*) ' Energy for upper state from exp. energy diff.'
      endif
      write(13,*) ' Energy (a.u.) infinite nucl. mass for upper state',
     1 eup,'H'

      write(13,*)
      write(*,*) ' Value of sms parameter S (a.u.) for lower state'
      read(*,*) gradlow
      write(13,*) ' Value of sms parameter S (a.u.) for lower state',
     1 gradlow,'H'
      
      write(*,*) ' Value of sms parameter S (a.u.) for upper state'
      read(*,*) gradup
      write(13,*) ' Value of sms parameter S (a.u.) for upper state',
     1 gradup,'H'

      write(13,*)
      write(*,*)  
     1    ' Electron density (a.u.^-3) at the nucleus for lower state'
      read(*,*) eldenslow
      write(13,*) 
     1    ' Electron density (a.u.^-3) at the nucleus for lower state',
     2       eldenslow
      
      write(*,*)  
     1    ' Electron density (a.u.^-3) at the nucleus for upper state'
      read(*,*) eldensup
      write(13,*) 
     1    ' Electron density (a.u.^-3) at the nucleus for upper state',
     2      eldensup

      write(*,*) ' Number of isotopes'
      read(*,*) niso 
      write(13,*)
      write(13,*) ' Number of isotopes',niso

      write(13,*)
      write(13,*) ' Isotope label, atomic isotope mass (u) and root mean square
     1 nuclear radius (fm)'
      do 10 i = 1,niso
         write(*,*) ' Isotope label (character*4)'
         read(*,'(a4)') isolabel(i)
         write(*,*) ' Atomic isotope mass (u)'
         read(*,*) isomass(i)
         write(*,*) ' Root meansquare nuclear radius sqrt(<r**2 >) (fm)'
         read(*,*) rms(i)
         write(13,*) 
     1   isolabel(i),'    ',isomass(i),'u','    ',rms(i),'fm'
10    continue 

Cmrg  write(*,*) ' Number of electrons' 
Cmrg  read(*,*) ne
Cmrg  write(13,*)
Cmrg  write(13,*) ' Number of electrons',ne 
Cmrg  Z*emass is the electron mass we need to substract from the
Cmrg  atomic mass to get the nuclear mass. Z is the number of protons,
Cmrg  ie. the number of electrons for the neutral atom.
Cmrg  M_N (A,Z) = M_A (A,Z) - Z*m_e + B_e(Z)  
C      ne = ifix(z)
      ne = idnint(z)
      write(*,*) 'Number of electrons',ne



*     Calculate the nuclear mass from atomic mass, that is subtract
*     the electron mass

      write(13,*)
      write(13,*) ' Nuclear mass for the isotopes'
      do 20 i = 1,niso
         isomass(i) = isomass(i) - ne*emass
         write(13,*) isolabel(i),'    ',isomass(i),'u'
         muovm(i) = isomass(i)/(isomass(i) + emass)
         write(13,*) '     mu/m ratio = ',muovm(i)
         rydm(i) = ry * muovm(i)
         write(13,*) '     Rydberg_M  = ',rydm(i)*0.01D0,' cm-1 '
         write(13,*) ' (2*)Rydberg_M  = ',
     :                                rydm(i)*0.02D0,' cm-1 '
         write(13,*)
         rme(i) = emass * muovm(i)
20    continue

      write(13,*)
      write(13,*) '============='
      write(13,*) ' Lower state'
      write(13,*) '============='


      write(13,*)
      write(13,*) ' Energy for infinite nuclear mass'

      write(13,*) '        ',elow,'H'

*     Lower energy for the different isotopes

      write(13,*)
      write(13,*) ' Energy corrected for normal mass effect'

      do 25 i = 1,niso
Cmrg     el(i) = elow*isomass(i)/(isomass(i) + emass) 
         el(i) = elow*muovm(i)
         eln(i) = el(i)
         write(13,*) isolabel(i),'    ', el(i),'H'
25    continue

      write(13,*)
      write(13,*)' Energy corrected for normal and specific mass effect'

      do 30 i = 1,niso
Cmrg     el(i) = el(i) + emass*gradlow/isomass(i)
         el(i) = el(i) + rme(i)*muovm(i)*gradlow/isomass(i)
         elns(i) = el(i)
         write(13,*) isolabel(i),'    ', el(i),'H'
30    continue        

      write(13,*)
      write(13,*)' Energy corrected for normal and specific mass effect'
     1          ,' and fieldshift'

      do 31 i = 1,niso
         rms(i) = rms(i)*1.d-15/a0
         el(i) = el(i) + 2.d0*pi*z*eldenslow*rms(i)*rms(i)/3.d0
         write(13,*) isolabel(i),'    ', el(i),'H'
31    continue        

      write(13,*)
      write(13,*) '============='
      write(13,*) ' Upper state'
      write(13,*) '============='

      write(13,*)
      write(13,*) ' Energy for infinite nuclear mass'

      write(13,*) '        ',eup,'H'

*     Upper energy for the different isotopes

      write(13,*)
      write(13,*) ' Energy corrected for normal mass effect'

      do 37 i = 1,niso
Cmrg     eu(i) = eup*isomass(i)/(isomass(i) + emass) 
         eu(i) = eup*muovm(i)
         eun(i) = eu(i)
         write(13,*) isolabel(i),'    ', eu(i),'H'
37    continue        

      write(13,*)
      write(13,*)' Energy corrected for normal and specific mass effect'

      do 40 i = 1,niso
Cmrg     eu(i) = eu(i) + emass*gradup/isomass(i)
         eu(i) = eu(i) + rme(i)*muovm(i)*gradup/isomass(i)
         euns(i) = eu(i)
         write(13,*) isolabel(i),'    ', eu(i),'H'
40    continue        

      write(13,*)
      write(13,*)' Energy corrected for normal and specific mass effect'
     1          ,' and fieldshift'

      do 41 i = 1,niso
         eu(i) = eu(i) + 2.d0*pi*z*eldensup*rms(i)*rms(i)/3.d0 
         write(13,*) isolabel(i),'    ', eu(i),'H'
41    continue        

      write(13,*)
      write(13,*) '===================='
      write(13,*) ' Energy differences'
      write(13,*) '===================='

      write(13,*)
      write(13,*) ' Energy difference between upper and lower state'
      write(13,*)
      write(13,*) ' Difference for infinite nuclear mass'
      write(13,*) 
     
*     Energy difference

      diff(1) = eup - elow
      write(13,*) '        ',diff(1),'H'

      write(13,*)
      write(13,*) '        ',diff(1)*ry*2.d0*0.01d0,'cm-1'

      write(13,*)
      write(13,*) ' Transition wavelength (Angstrom)'

      wlength(1) = 1.d10/(2.d0*diff(1)*ry)
      write(13,*) '        ', wlength(1),'A'

      write(13,*)
      write(13,*) ' Differences with normal mass correction'

*     Energy difference

      do 50 i = 1,niso
         diff(i) = eun(i) - eln(i)
         write(13,*) isolabel(i),'    ',diff(i),'H'
50    continue

      write(13,*)
      do 55 i = 1,niso
         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
55    continue

      write(13,*)
      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
      dif1 = diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0
      

      write(13,*)
      write(13,*) ' Transition wavelength (Angstrom)'

      do 60 i = 1,niso
         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
         write(13,*) isolabel(i),'    ', wlength(i),'A'
60    continue

      write(13,*)
      write(13,*) ' Wavelength difference between the isotopes'

      do 70 i = 2,niso
         dwlength(i) = dabs(wlength(1) - wlength(i))
         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
70    continue

      write(13,*)
      write(13,*) 
     :' Differences with normal and specific mass correction'

*     Energy difference

      do 150 i = 1,niso
         diff(i) = euns(i) - elns(i)
         write(13,*) isolabel(i),'    ',diff(i),'H'
150    continue

      write(13,*)
      do 155 i = 1,niso
         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
155    continue

      write(13,*)
      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
      dif2 = diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0

      dif21 = dif2 - dif1
      write(13,*)
      write(13,*) ' Transition SMS (cm-1) = ', dif21
      dif21 = dif21 * 2.99792458D+04
      write(13,*) ' Transition SMS (MHz)  = ', dif21


      write(13,*)
      write(13,*) ' Transition wavelength (Angstrom)'

      do 160 i = 1,niso
         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
         write(13,*) isolabel(i),'    ', wlength(i),'A'
160    continue

      write(13,*)
      write(13,*) ' Wavelength difference between the isotopes'

      do 170 i = 2,niso
         dwlength(i) = dabs(wlength(1) - wlength(i))
         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
170    continue

      write(13,*)
      write(13,*)
      write(13,*) ' Differences with normal and specific mass',
     1            ' and field shift correction'

*     Energy difference

      do 180 i = 1,niso
         diff(i) = eu(i) - el(i)
         write(13,*) isolabel(i),'    ',diff(i),'H'
180    continue

      write(13,*)
      do 185 i = 1,niso
         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
185    continue

      write(13,*)
      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
      write(13,*)

      write(13,*) ' Transition wavelength (Angstrom)'

      do 190 i = 1,niso
         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
         write(13,*) isolabel(i),'    ', wlength(i),'A'
190    continue

      write(13,*)
      write(13,*) ' Wavelength difference between the isotopes'

      do 200 i = 2,niso
         dwlength(i) = dabs(wlength(1) - wlength(i))
         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
200    continue

      end
