##
## Makefile for testing module
##

DEBUG 		= Yes

MKDIR		= mkdir -p

PARENT_BIN	= $(shell pwd)/../bin

SRC_DIR		= $(shell pwd)/src

BIN_DIR		= $(shell pwd)/bin

LINKERS		= -I ${PARENT_BIN}
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

TARGET		= $(BIN_DIR)/calculate_error.out

MODULES		= ${PARENT_BIN}/kinds.mod \
			  ${PARENT_BIN}/constants.mod \
			  ${PARENT_BIN}/essentials.mod
			 
OBJECTS		= ${PARENT_BIN}/kinds.o \
			  ${PARENT_BIN}/constants.o \
			  ${PARENT_BIN}/essentials.o \
			  $(BIN_DIR)/calculate_error.o

all : $(BIN_DIR) $(TARGET)


#
# The "clean" target just removes all files (modules, object code)
# produced during the build process:
#
.PHONY: clean
clean::
		$(RM) -rf $(TARGET) ${BIN_DIR}/calculate_error.o

#
# Rule to create directories if not present
#

$(BIN_DIR) :
	$(MKDIR) $(BIN_DIR)

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
${BIN_DIR}/%.o: ${SRC_DIR}/%.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -o $@ -c $< $(LINKERS)
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -o $@ -c $< $(LINKERS)
endif

#
# Pattern rule that produces a Fortran module from source
# code:
#
${BIN_DIR}/%.mod: ${SRC_DIR}/%.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -c $< $(LINKERS)
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -c $< $(LINKERS)
endif

