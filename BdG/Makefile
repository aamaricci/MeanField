##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =.

##$ SET THE TARGET DIRECTORY WHERE TO PUT THE EXECUTABLE (default if $HOME/.bin in the PATH)
DIREXE=$(HOME)/.bin

##$ CHOOSE THE DRIVER CODE:
EXE=bdg_disorder

##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)


ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif

ifeq ($(PLAT),gnu)
NFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
FFLAG = -O3 -ffast-math -march=native -funroll-all-loops -ffree-line-length-none
FPPFLAG =-cpp
endif



##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files

.SUFFIXES: .f90



###############################################################################
#                        HERE STARTS THE REAL WORK
###############################################################################

OBJS= 

all: all version compile completion
mpi: mpi version compile completion
debug: debug version compile completion
debug_mpi:debug_mpi version compile completion

all: FPPFLAG+=-D_

mpi: FPPFLAG+=-D_MPI

debug: FFLAG=$(DFLAG)
debug: FPPFLAG+=-D_

debug_mpi: FFLAG=${DFLAG} 
debug_mpi: FPPFLAG+=-D_MPI


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FPPFLAG) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(INCARGS) $(LIBARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) -c $<

completion:
	scifor_completion.sh $(DIR)/$(EXE).f90

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)


