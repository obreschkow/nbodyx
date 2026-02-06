# Call as make [mode=standard,dev] [kind=8,16]

# Compiler settings:
# mode = compilation mode; allowed modes are 'standard' and 'dev'. The latter is a developer mode that runs slower with more checks.
# kind = 8 or 16; number of bytes used to internally represent particle positions; 16 is only needed for very high accuracy
#        computations, coupled with a high-order integrator and small time-steps. On standard 64-bit systems, simulations run with
#        16 bytes take about twice the computation time.
#        Note: make clean is required when switching between compiler settings (mode and kind).

ifndef mode
   mode = standard
endif

ifndef kind
   kind = 8
endif

# standard compiler flags (depend on the "mode" option)
ifeq ($(mode),standard)
   CFLAGS = -Ofast -fdefault-real-8
else ifeq ($(mode),dev)
   CFLAGS = -O0 -g -fbounds-check -fwhole-file -ffpe-trap=invalid,zero,overflow -Wall -Wunused -Wuninitialized -Wsurprising -Wconversion
else
   $(info ERROR unknown mode: '${mode}')
stop
endif

# concatenate flags
FCFLAGS =  $(CFLAGS)

# Compiler
FC = gfortran

# user info
$(info Compilation options:)
$(info + Compiling mode = '${mode}'.)

# List of executables to be built within the package
PROGRAMS = nbodyx

# "make" builds all
all: $(PROGRAMS)

nbodyx.o:   shared_module_core.o \
				shared_module_arguments.o \
				shared_module_parameters.o \
				shared_module_constants.o \
				module_particle_$(kind).o \
				module_global.o \
				module_io.o \
				module_acceleration.o \
				module_integration.o
nbodyx: 	   shared_module_core.o \
				shared_module_arguments.o \
				shared_module_parameters.o \
				shared_module_constants.o \
				module_particle_$(kind).o \
				module_global.o \
				module_io.o \
				module_acceleration.o \
				module_integration.o

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
	
%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) -c $<
	
%.o: %.f08
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	rm -f *~ $(PROGRAMS)
	rm -f fort.*
	rm -rf output