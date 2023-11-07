################################################################################
#  Makefile for struture grid package

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC    := $(GNUHOME)/bin/gcc
CXX    := $(GNUHOME)/bin/g++

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./src/ $(CFLAGS)

#- O3
CFLAGS := -O3 -std=c99 $(CFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
#LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

grid_post_proc_3d: cJSON.o lib_mem.o \
				io_funcs.o gd_t.o par_t.o \
				quality_check.o algebra.o main.o
	$(CC) -o $@ $^ $(LDFLAGS)

cJSON.o: lib/cJSON.c
	${CC} -c -o $@ $(CFLAGS) $<
lib_mem.o: lib/lib_mem.c
	${CC} -c -o $@ $(CFLAGS) $<
par_t.o: src/par_t.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_t.o: src/gd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
io_funcs.o: src/io_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
algebra.o: src/algebra.c
	${CC} -c -o $@ $(CFLAGS) $<
quality_check.o:src/quality_check.c
	${CC} -c -o $@ $(CFLAGS) $<

main.o: src/main.c
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f grid_post_proc_3d

cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
