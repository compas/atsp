
################################################################################
# Makefile for all applications. 
# .For each subdirectory the default target "install" is recursively built 
# ."make clean" will descend into each directory and remove *.o files.
################################################################################


TARGETS= install 
SUBDIRS= bp_eiv bp_ang bp_mat bpci hfs iso lsreduce mchf nonh trans biotr \
         mpi_bp_ang mpi_bp_eiv mpi_bp_mat mpi_mchf \
         mpi_nonh mpi_bpci


$(TARGETS): $(SUBDIRS)
install:
	@for i in $(SUBDIRS) ; \
	do \
                cd $$i ; \
                echo "....Entering: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
                echo "....Leaving: $$i" ; echo; echo;\
	done

backup:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
		$(MAKE) $@ ; \
		cd .. ; \
	done

clean:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
	        $(MAKE) $@ ; \
	        cd .. ; \
	done

