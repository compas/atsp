# makefile for Applications

TARGETS= install #debug
SUBDIRS= libang libcom libdvdson librad libnet mpi

$(TARGETS) : $(SUBDIRS)
install:
	@for i in $(SUBDIRS) ; \
	do \
                cd $$i ; \
                echo "....Entering: $$i" ; \
                $(MAKE) $@ ; \
                cd .. ; \
                echo "....Leaving: $$i" ; echo; echo;\
	done

clean:
	@for i in $(SUBDIRS); \
	do \
	cd $$i ; \
	echo "Current directory: $$i" ; \
		$(MAKE) $@ ; \
	        cd .. ; \
	done

