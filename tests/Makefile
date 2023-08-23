include ../Make.inc

SP3MM_DIR=..
SP3MM_INC_DIR=$(SP3MM_DIR)/include
SP3MM_MOD_DIR=$(SP3MM_DIR)/modules
SP3MM_LIB_DIR=$(SP3MM_DIR)/lib
SP3MM_LIBS=-L$(SP3MM_LIB_DIR) -lsp3mm
EXEC_DIR=runs
RESULT_DIR=results

FINCLUDES=-I. -I$(SP3MM_INC_DIR) -I$(SP3MM_MOD_DIR) $(PSBLAS_INCLUDES)
LINKOPT=-fopenmp

TESTS=	sp3mm_test_serial \
		sp3mm_test_serial_reversed \
		sp3mm_test_omp_gustavson \
		sp3mm_test_omp_gustavson_reversed \
		sp3mm_test_omp_gustavson_1d \
		sp3mm_test_omp_gustavson_1d_reversed \
		sp3mm_test_serial_rb_tree \
		sp3mm_test_serial_rb_tree_reversed \
		sp3mm_test_omp_rb_tree \
		sp3mm_test_omp_rb_tree_reversed \
		sp3mm_test_omp_two_pass \
		sp3mm_test_omp_two_pass_reversed

all: exec_dir result_dir $(TESTS)

exec_dir:
	mkdir -p $(EXEC_DIR)

result_dir:
	mkdir -p $(EXEC_DIR)/$(RESULT_DIR)

sp3mm_test_%: sp3mm_test_mod.o sp3mm_test_%.o
	$(FLINK) $(LINKOPT) $^ -o $@ $(SP3MM_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv $@ $(EXEC_DIR)

clean:
	/bin/rm -f *.o
	/bin/rm -f runs/sp3mm_test*
	/bin/rm -f runs/test_settings.conf
