##
## Makefile for testing module
##

DEBUG 		= Yes

BIN_DIR		= $(shell pwd)/..

LINKERS		= -I ${BIN_DIR}
#
# Fortran Compiler (FC) to be used:
#
FC		= ifort

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

TARGET		= calculate_error.out

MODULES		= ${BIN_DIR}/kinds.mod \
			  ${BIN_DIR}/constants.mod \
			  ${BIN_DIR}/essentials.mod
			 
OBJECTS		= ${BIN_DIR}/kinds.o \
			  ${BIN_DIR}/constants.o \
			  ${BIN_DIR}/essentials.o \
			  calculate_error.o

all : $(TARGET)


#
# The "clean" target just removes all files (modules, object code)
# produced during the build process:
#
.PHONY: clean
clean::
		$(RM) -rf $(TARGET) calculate_error.o

#
# Rule to produce the target executable:
#
	
$(TARGET) : $(MODULES) $(OBJECTS)
ifeq ($(DEBUG), Yes)
	$(FC) $(FFLAGS_FOR_DEBUG) -o $@ $(LINKERS) $(OBJECTS) $(LDFLAGS) $(LIBS)
else
	$(FC) $(FFLAGS_FOR_RUN) -o $@ $(LINKERS) $(OBJECTS) $(LDFLAGS) $(LIBS)
endif

#
# Pattern rule that produces object code from source code:
#
%.o: %.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -o $@ -c $< $(LINKERS)
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -o $@ -c $< $(LINKERS)
endif

#
# Pattern rule that produces a Fortran module from source
# code:
#
%.mod: %.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -c $< $(LINKERS)
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -c $< $(LINKERS)
endif

