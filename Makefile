##
## Makefile for FDPML project
##

DEBUG 		= Yes

MKDIR		= mkdir -p

SRC_DIR		= $(shell pwd)/src

BIN_DIR		= $(shell pwd)/bin

#
# List of Fortran modules that should be generated from source files
# (with .f90 suffix) in this directory
#
MODULES		= ${BIN_DIR}/kinds.mod \
			${BIN_DIR}/init.mod \
			${BIN_DIR}/constants.mod \
			${BIN_DIR}/mp_module.mod \
			${BIN_DIR}/essentials.mod \
			${BIN_DIR}/io_module.mod \
			${BIN_DIR}/rigid.mod \
			${BIN_DIR}/ws.mod \
			${BIN_DIR}/sum_rule.mod \
			${BIN_DIR}/matrix_inversion.mod \
			${BIN_DIR}/monkhorstpack.mod \
			${BIN_DIR}/dispersion.mod \
			${BIN_DIR}/gnufor2.mod \
			${BIN_DIR}/COO_routines.mod \
			${BIN_DIR}/preprocessing_module.mod \
			${BIN_DIR}/matgen_module.mod \
			${BIN_DIR}/bicg.mod 


#
# Sources that go into the executable (not modules):
#
SOURCES		= FDPML.f90

#
# Generate list of object files for the SOURCES listed
# above:
#
OBJECTS		= ${BIN_DIR}/kinds.o \
			${BIN_DIR}/init.o \
			${BIN_DIR}/constants.o \
			${BIN_DIR}/mp_module.o \
			${BIN_DIR}/essentials.o \
			${BIN_DIR}/io_module.o \
			${BIN_DIR}/rigid.o \
			${BIN_DIR}/ws.o \
			${BIN_DIR}/sum_rule.o \
			${BIN_DIR}/matrix_inversion.o \
			${BIN_DIR}/monkhorstpack.o \
			${BIN_DIR}/dispersion.o \
			${BIN_DIR}/gnufor2.o \
			${BIN_DIR}/COO_routines.o \
			${BIN_DIR}/preprocessing_module.o \
			${BIN_DIR}/matgen_module.o \
			${BIN_DIR}/bicg.o \
			${BIN_DIR}/FDPML.o

#
# What should the executable produced be named?
#
TARGET		= ${BIN_DIR}/FDPML.out

##
####
####
##

#
# Fortran Compiler (FC) to be used:
#
FC		= mpifort

#
# Fortran FLAGS (FFLAGS) to be used (not preprocessor definitions,
# just flags that affect compilation dialect, optimizations, etc.):
#
# compiler flags for final executable
FFLAGS_FOR_RUN		= -mkl 

# compiler flags for execulatable with call-back and traceback capabilities 
#(add other flags that you might require for debugging)
FFLAGS_FOR_DEBUG	= -mkl -CB -traceback

#
# Fortran Pre-Processor FLAGS (FPPFLAGS) that affect compilation
# (e.g. "-I../common" or "-DPI=3.14159"):
#
FPPFLAGS	= 

#
# LD FLAGS (LDFLAGS) are flags used during the linking stage that
# are specific to the linker:
#
LDFLAGS		= 

#
# Any libraries needed by the linker (e.g. "-lm"):
#
LIBS		= -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64

##
####
####
##

#
# The default target if nothing is supplied with the "make"
# command:
#
default: $(TARGET)

#
# Run "make all" command to generate target file
#
all : 	$(BIN_DIR) $(TARGET)
		+$(MAKE) -C tests

#
# The "clean" target just removes all files (modules, object code)
# produced during the build process:
#
.PHONY: clean
clean::
		$(RM) -rf $(TARGET) $(MODULES) $(OBJECTS)
		+$(MAKE) -C tests clean
		
#
# Rule to create directories if not present
#

$(BIN_DIR) :
	$(MKDIR) $(BIN_DIR)


#
# Rule to produce the target executable:
#
$(TARGET): $(OBJECTS)
ifeq ($(DEBUG), Yes)
	$(FC) $(FFLAGS_FOR_DEBUG) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)
else
	$(FC) $(FFLAGS_FOR_RUN) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)
endif
#
# Pattern rule that produces object code from source code:
#
${BIN_DIR}/%.o: ${SRC_DIR}/%.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -module $(BIN_DIR) -c -o $@ $<
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -module $(BIN_DIR) -c -o $@ -c $<
endif

#
# Pattern rule that produces a Fortran module from source
# code:
#
#${BIN_DIR}/%.mod: ${SRC_DIR}/%.f90
#ifeq ($(DEBUG), Yes)
#	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -c -module $(BIN_DIR) -o $< 
#else
#	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -c -module $(BIN_DIR) -o $< 
#endif


