# ----------------------------------------------------------------------------------------
# EDIT THESE
COMPILER := intel
FC := mpif90
SRCDIR := source
BUILDDIR := build
BINDIR := bin
EXECUTABLE := $(BINDIR)/lesmpi.exe

# ----------------------------------------------------------------------------------------
# DON'T TOUCH THESE
$(shell mkdir -p $(BUILDDIR)) # make the build directory if it doesn't exist

# `shell find` all directories, `filter-out` the SRCDIR itself, `sort` alphabetically
SRC_DIRS := $(sort $(filter-out $(SRCDIR), $(shell find $(SRCDIR) -type d)))

# `shell find` gets all filenames, including directory paths, then `notdir` strips off directory paths
MOD_SRCS := inputs.f90 pars.f90 con_data.f90 con_stats.f90 fftwk.f90 fields.f90 tracerbc.f90 reaction.f90
F90_SRCS := $(notdir $(shell find $(SRCDIR) -name '*.f90'))
F90_SRCS := $(filter-out $(MOD_SRCS), $(F90_SRCS)) # strip out modules from other f90 files
F_SRCS := $(notdir $(shell find $(SRCDIR) -name '*.f'))

# Generate corresponding object file paths in the build directory
MOD_OBJS := $(MOD_SRCS:%.f90=$(BUILDDIR)/%.mod)
F90_OBJS := $(F90_SRCS:%.f90=$(BUILDDIR)/%.o)
F_OBJS := $(F_SRCS:%.f=$(BUILDDIR)/%.o)

vpath %.f90 $(SRC_DIRS) # vpath tells make where to look for files without explicit paths
vpath %.f $(SRC_DIRS)

# include all the compiler-specific options (necessary)
mkfile := make.$(strip $(COMPILER))
$(info Loading compiler options from $(mkfile))
include $(mkfile)

# ----------------------------------------------------------------------------------------
# EDIT THESE
# OPTIONS here are defined in the make.$(COMPILER) files

# As the first recipe in this list, running just `make` aliases to `make syntax`
# Rule for just checking syntax bugs, same kind of info as if vscode were linting your code for you
.PHONY: syntax
syntax: OPTIONS=$(SYNTAX)
syntax: $(MOD_OBJS) $(F90_OBJS)

.PHONY: debug
debug: OPTIONS=$(DBG1)
debug: $(EXECUTABLE)

.PHONY: profile
profile: OPTIONS=$(OPT2)
profile: $(EXECUTABLE)

.PHONY: fast
fast: OPTIONS=$(OPT4)
fast: $(EXECUTABLE)

# aliases the common `make all` command to `make fast`
all: fast

# a command to reformat files
RAW_SRCS := $(shell find $(SRCDIR) -name '*.f90') # includes directory paths
.PHONY: format
format:
	fprettify --indent 4 -w 4 -l 300 --whitespace-intrinsics true\
	 --enable-decl --enable-replacements --c-relations $(RAW_SRCS)

# This removes the build directory and executable in teh test directory
.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(EXECUTABLE) $(EXECUTABLE).dSYM

# this removes everything in the .build/ directory AND any accidental outputs made in the project directory
.PHONY: realclean
realclean: clean
	rm -rf *.o *.mod *.dSYM

# ----------------------------------------------------------------------------------------
# DON'T THOUCH THESE
# MODFLAGS, F90FLAGS, FFLAGS here are defined in the make.$(COMPILER) files

# Rule to link all .mod and .o object binaries into a single program
$(EXECUTABLE): $(MOD_OBJS) $(F_OBJS) $(F90_OBJS)
	@mkdir -p $(BINDIR)
	$(FC) $(F90FLAGS) $(OPTIONS) $(BUILDDIR)/*.o -o $@ $(LDFLAGS)

# Rule to compile module .f90 files into .mod files
$(BUILDDIR)/%.mod: %.f90
	$(FC) $(F90FLAGS) $(OPTIONS) -c $< -o $(BUILDDIR)/$*.o

# Rule to compile .f files into .o object binaries
$(BUILDDIR)/%.o: %.f
	$(FC) $(FFLAGS) $(OPTIONS) -c $< -o $@

# Rule to compile non-module .f90 files into .o object binaries
$(BUILDDIR)/%.o: %.f90
	$(FC) $(F90FLAGS) $(OPTIONS) -c $< -o $@
