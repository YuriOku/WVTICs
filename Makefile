SHELL = /bin/bash

## OPTIONS  ##
# always recommended
OPT     += -DOUTPUT_DIAGNOSTICS     # Output extra diagnostics to file
OPT     += -DGADGET2_OUTPUT_ORDER

# need for 2D
#OPT     += -DTWO_DIM               # 2D instead of 3D -> set z component to 0

# not needed
#OPT     += -DSAVE_WVT_STEPS         # write IC file for every WVT step
#OPT     += -DSPH_CUBIC_SPLINE      # for use with Gadget2
#OPT     += -DSPH_WC2               # wendland c2 kernel
#OPT     += -DEAT_PNG                 # Eat density profile from a png file
#OPT     += -DPEANO_SAMPLING      # use peano curve based sampling to improve initial random positions
#OPT     += -DBRUTE_FORCE_NGB               # Use a brute force neighbour finder instead of the tree based one
#OPT     += -DUSE_APM               # (experimental) Use artificial pressure method (Rosswog 2020)

# not recommended
#OPT     += -DREJECTION_SAMPLING      # use von Neumann rejection sampling to improve initial random positions
#OPT     += -DREGULARIZE             # regularize particle distribution in a similar way with Gaburov & Nitadori (2011)

ifndef SYSTYPE
    SYSTYPE := $(shell hostname)
endif

CC       = icc
OPTIMIZE = -Wall -O3 -xHost
GSL_INCL = -I/home/yuri/local/include
GSL_LIBS = -L/home/yuri/local/lib
PNG_LIBS = -L/usr/local/Cellar/libpng/1.6.21/lib -lpng -lz
PNG_INCL = -I/usr/local/Cellar/libpng/1.6.21/include


## TARGET ##

EXEC = WVTICs

## FILES ##

SRCDIR    = src/

SRCFILES := ${shell find $(SRCDIR) -name \*.c -print} # all .c files in SRCDIR
OBJFILES = $(SRCFILES:.c=.o)

INCLFILES := ${shell find src -name \*.h -print} # all .h files in SRCDIR
INCLFILES += Makefile

CFLAGS  = -std=c99 -qopenmp $(OPTIMIZE) $(OPT) $(GSL_INCL) $(PNG_INCL)

LINK    = $(GSL_LIBS) -lm -lgsl -lgslcblas $(PNG_LIBS)

## RULES ##

%.o : %.c
	@echo [CC] $@
	@$(CC) $(CFLAGS)  -o $@ -c $<
#	$(CC) $(CFLAGS)  -o $@ -c $<

$(EXEC) : $(OBJFILES)
	$(CC) $(CFLAGS) $(OBJFILES) $(LINK) -o $(EXEC)
#	@ctags -w $(SRCFILES) $(INCLFILES)

$(OBJFILES) : $(INCLFILES) $(SRCFILES)

clean :
	rm $(OBJFILES) $(EXEC)
