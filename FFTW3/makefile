include ../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src
INCDIR = inc
LIBDIR = lib

EXE = prg

ifeq ($(FORT),x86_64-w64-mingw32-gfortran-win32)
	LIBDIR = libwin64
	EXE = prg.exe
endif

ALL_LIBS = $(addprefix $(LIBDIR)/, *.a)

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ       = $(SRC:.$(EXT)=.o)
ALL_OBJS  = $(addprefix $(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(utils_dir)/$(OBJDIR)/, mod_*.o)

ALL_MODS  = -I$(MODDIR) -I$(INCDIR)
ALL_MODS += -I$(utils_dir)/$(MODDIR)

CFLAGS  = $(ALL_MODS) -fopenmp -fPIC
CFLAGS += -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none

LFLAGS  = $(ALL_LIBS)
LFLAGS += -lm -lgomp -lpthread
# be careful, -lpthread must be at the end

ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -$(OPTC)
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

mod_fftw3.o :
prg.o : mod_fftw3.o

#--------------------------------------------------
.PHONY: clean debug gprof all last

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

