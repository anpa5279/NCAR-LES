# MODIFY THESE AS YOU NEED/LIKE

COMMON := -fdefault-real-8 -fdefault-double-8 -fPIC -pipe -J./build -std=legacy
FFLAGS := $(COMMON) -ffixed-line-length-132
F90FLAGS := $(COMMON) -ffree-line-length-none

LDFLAGS := -lm

SYNTAX := -fsyntax-only -Wall # -Wfatal-errors
REFACTORING := -fsyntax-only -fimplicit-none -Wall -Wextra -Wfatal-errors -std=f03

# Og allows some O1 optimizations while making debugging cleaner/clearer than O0
# also, plain `-g` is equivalent to `-g2`
DBG1 := -Og -g1 -ffpe-trap=invalid,zero,overflow
DBG2 := -Og -g2 -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all,no-array-temps -fmax-errors=5
DBG3 := -O0 -g2 -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all -fmax-errors=5

# CODE PROFILING
warn := -Warray-temporaries -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure\
		-Wrealloc-lhs-all -Wfrontend-loop-interchange
WRN1 := -g2 -O1 $(warn)
WRN2 := -g2 -O2 $(warn)

OPT1 := -g2 -O2
OPT2 := -g2 -O2 -ffast-math -fno-protect-parens

# HARDCORE OPTIMIZATION
OPT3 := -O3
OPT4 := -O3 -ffast-math -fno-protect-parens
