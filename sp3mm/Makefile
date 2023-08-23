include ../Make.inc

LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
HERE=.

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)

MODOBJS= sp3mm_idx_tree_mod.o spmm_mod.o

OBJS=$(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libsp3mm.a

objs: $(OBJS)
	/bin/cp -p *$(.mod) $(MODDIR)

lib: objs
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)

$(MODOBJS): $(PSBLAS_MODDIR)/$(PSBBASEMODNAME)$(.mod)

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)