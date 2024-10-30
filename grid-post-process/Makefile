################################################################################
#  Makefile for struture grid package

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC    :=  $(GNUHOME)/bin/gcc

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./src/ $(CFLAGS)

#- O3
CFLAGS := -O3 -std=c99 $(CFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
#LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

skeldirs := obj
DIR_OBJ  := ./obj
#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

OBJS := cJSON.o lib_mem.o lib_math.o\
			  io_funcs.o gd_t.o par_t.o \
				quality_check.o algebra.o main.o

OBJS := $(addprefix $(DIR_OBJ)/,$(OBJS))

vpath  %.c .

all: skel main
skel:
	@mkdir -p $(skeldirs)

main: $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

$(DIR_OBJ)/%.o : src/%.c
	${CC} $(CFLAGS) -c $^ -o $@ 
$(DIR_OBJ)/%.o : lib/%.c
	${CC} $(CFLAGS) -c $^ -o $@ 

cleanexe:
	rm -f main

cleanobj:
	rm -rf $(DIR_OBJ)
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
