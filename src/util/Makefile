.SUFFIXES: .f90 .mod

# Flags to be added for compilation
NOFFLAGS = 

OBJU = \
comp.o        plotw.o       ondens.o     hf.o         \
printw.o      w_format.o    gencl.o       levels.o      \
relabel.o     w_unformat.o  lines.o       lsgen.o\
select.o      csfexcitation.o  isosupport2.o

SRCU = \
comp.f        plotw.f       condens.f     hf.f         \
printw.f      w_format.f    gencl.f       levels.f      \
relabel.f     w_unformat.f  lines.f       lsgen.f   \
select.f      csfexcitation.f90  isosupport2.f


LIBDIR = $(ATSP)/lib
LIBUSE = -lcom -lrad -lang -lcom

install : hf  comp condens  levels lines \
	  plotw printw relabel w_format   \
          w_unformat gencl select lsgen  csfexcitation isosupport2
	
	cp save $(ATSP)/bin
	cd tbl; make; cd ..
	cd LS_trends; make; cd ..

hf: hf.o
	$(FC)  -o hf $(FC_FLAGS) hf.o
	mv hf $(ATSP)/bin

lsgen : lsgen.o
	$(FC) -o lsgen $(NOFFLAGS) -L$(LIBDIR) lsgen.f $(LIBUSE)
	mv lsgen $(ATSP)/bin

comp: comp.o
	$(FC)  -o comp $(FC_FLAGS)  comp.o 
	mv comp $(ATSP)/bin

condens : condens.o 
	$(FC) -o condens $(FC_FLAGS) -L$(LIBDIR) condens.o  $(LIBUSE)
	mv  condens $(ATSP)/bin

gencl : gencl.o
	${FC} -o gencl $(FC_FLAGS) -L$(LIBDIR) gencl.f $(LIBUSE)
	mv gencl $(ATSP)/bin


levels : levels.o
	$(FC) -o levels $(FC_FLAGS) -L$(LIBDIR) levels.o $(LIBUSE)
	mv levels $(ATSP)/bin

lines: lines.o
	$(FC)  -o lines $(FC_FLAGS) lines.o
	mv lines $(ATSP)/bin

printw: printw.o
	$(FC)  -o  printw $(FC_FLAGS) printw.o
	mv  printw $(ATSP)/bin

plotw : plotw.o
	$(FC)  -o  plotw $(FC_FLAGS) plotw.o
	mv  plotw $(ATSP)/bin

select: select.o 
	$(FC) -o select $(FC_FLAGS) -L$(LIBDIR) select.o  $(LIBUSE)
	mv  select $(ATSP)/bin

csfexcitation: csfexcitation.o
	$(FC) -o csfexcitation $(FC_FLAGS) csfexcitation.o 
	mv csfexcitation $(ATSP)/bin

relabel : relabel.o
	$(FC) -o relabel $(FC_FLAGS)  relabel.o
	mv relabel $(ATSP)/bin

w_format : w_format.o
	$(FC) -o w_format $(FC_FLAGS)  w_format.o
	mv w_format $(ATSP)/bin

w_unformat : w_unformat.o
	$(FC) -o w_unformat $(FC_FLAGS)  w_unformat.o
	mv w_unformat $(ATSP)/bin

isosupport2 : isosupport2.o
	$(FC) -o TIS $(FC_FLAGS)  isosupport2.o
	mv TIS $(ATSP)/bin

clean:
	-rm -f *.o

%.o : %.f
	$(FC) -c $(FC_FLAGS) $< -o $@

%.o : %.f90
	$(FC) -c $(FC_FLAGS) $< -o $@

%.o : %.cc
	$(CPP) -c $(CPP_FLAGS) $< -o $@
