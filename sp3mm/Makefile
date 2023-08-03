include ../Make.inc

LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
HERE=.

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)

MODOBJS=sp3mm_const_mod.o sp3mm_config_mod.o \
		sp3mm_utils_mod.o sp3mm_acc_mod.o sp3mm_idx_tree_mod.o \
		spmm_mod.o sp3mm_base_mod.o

OBJS=$(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libsp3mm.a

all: objs impld

objs: $(OBJS)
	/bin/cp -p *$(.mod) $(MODDIR)

impld: objs
	$(MAKE) -C impl

lib: objs impld
	$(MAKE) -C impl lib
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)

$(MODOBJS): $(PSBLAS_MODDIR)/$(PSBBASEMODNAME)$(.mod)

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: implclean
	/bin/rm -f $(OBJS) *$(.mod)

implclean:
	$(MAKE) -C impl clean