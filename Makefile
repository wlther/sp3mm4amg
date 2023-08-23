include Make.inc

all: objs lib
	@echo ================================
	@echo Compilation of sp3mm succesfull
	@echo ================================

objs: sp3mm_objs

lib: libdir objs
	$(MAKE) -C sp3mm lib

libdir:
	mkdir -p include
	mkdir -p lib
	mkdir -p modules

sp3mm_objs:
	$(MAKE) -C sp3mm objs

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd modules; /bin/rm -f *.a *$(.mod) *$(.fh))

veryclean: cleanlib
	$(MAKE) -C sp3mm veryclean
	$(MAKE) -C tests veryclean

clean:
	$(MAKE) -C sp3mm clean