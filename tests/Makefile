include ../Make.inc

SP3MM_DIR=..
SP3MM_INC_DIR=$(SP3MM_DIR)/include
SP3MM_MOD_DIR=$(SP3MM_DIR)/modules
SP3MM_LIB_DIR=$(SP3MM_DIR)/lib
SP3MM_LIBS=-L$(SP3MM_LIB_DIR) -lsp3mm
EXEC_DIR=runs

FINCLUDES=-I. -I$(SP3MM_INC_DIR) -I$(SP3MM_MOD_DIR) $(PSBLAS_INCLUDES) -I/opt/papi-7.0.1/include
LINKOPT=-fopenmp

all: exec_dir sp3mm_test

exec_dir:
	mkdir -p $(EXEC_DIR)

sp3mm_test: sp3mm_test_mod.o sp3mm_test.o
	$(FLINK) $(LINKOPT) $^ -o $@ $(SP3MM_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv $@ $(EXEC_DIR)

clean:
	/bin/rm -f *.o
	/bin/rm -f runs/sp3mm_test

.PHONY: sp3mm_test