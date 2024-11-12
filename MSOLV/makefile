include ../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
INCDIR = inc

GENLIB  =
GENLIB += $(LIBDIR)/libopenblas.a

MUMLIB  =
MUMLIB += $(LIBDIR)/libmpiseq.a
MUMLIB += $(LIBDIR)/libdmumps.a
MUMLIB += $(LIBDIR)/libmumps_common.a
MUMLIB += $(LIBDIR)/libesmumps.a
MUMLIB += $(LIBDIR)/libpord.a

MUMLIB += $(LIBDIR)/libmetis.a
MUMLIB += $(LIBDIR)/libscotcherr.a
MUMLIB += $(LIBDIR)/libscotch.a

UMFLIB  =
UMFLIB += $(LIBDIR)/libumfpack.a
UMFLIB += $(LIBDIR)/libamd.a
UMFLIB += $(LIBDIR)/libsuitesparseconfig.a
UMFLIB += $(LIBDIR)/libcholmod.a
UMFLIB += $(LIBDIR)/libcolamd.a
UMFLIB += $(LIBDIR)/libccolamd.a
UMFLIB += $(LIBDIR)/libcamd.a

SUPLIB  =
SUPLIB += $(LIBDIR)/libsuperlu.a

M48LIB  =
ifeq ($(MA48_LIB),1)
   M48LIB += $(LIBDIR)/libhsl_ma48.a
endif

SMALIB = -Wl,--start-group $(GENLIB) $(M48LIB) $(MUMLIB) $(SUPLIB) $(UMFLIB) $(FFTLIB) -Wl,--end-group

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = prg

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ      = $(SRC:.$(EXT)=.o)
ALL_OBJS = $(addprefix $(OBJDIR)/, *.o) $(addprefix $(utils_dir)/$(OBJDIR)/, *.o) $(addprefix $(qsort_dir)/$(OBJDIR)/, *.o) 


CFLAGS  = -I$(MODDIR) -I$(INCDIR) -I$(qsort_dir)/$(MODDIR) -I$(utils_dir)/$(MODDIR)
CFLAGS += -cpp -DWITH_MA48=$(MA48_LIB) -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -fopenmp -fPIC

LFLAGS  = $(SMALIB)
LFLAGS += -lpthread -lm -lgomp


ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -ffpe-summary=invalid,zero,overflow,underflow -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
else
	CFLAGS += -O3 -ffpe-summary=none
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

%.o:	%.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	@find . -maxdepth 1 -name '*.mod*' -type f -print0 | xargs -0r mv -t ./$(MODDIR)

$(EXE):	$(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)
	@rm $(OBJDIR)/prg.o

mod_gen_par.o :

umfpack.o :
superlu.o :
dmumps_struc.o :
hsl_common90.o :
hsl_ddeps90.o :
hsl_ma48d.o : hsl_common90.o hsl_ddeps90.o

mod_solver.o : umfpack.o hsl_ma48d.o superlu.o dmumps_struc.o mod_gen_par.o

prg.o : mod_solver.o

#--------------------------------------------------
.PHONY: clean debug gprof all

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f $(EXE)

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE"

all:
	(make clean ; make)

