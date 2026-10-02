# SpheroDSMC: zero-dimensional DSMC of inelastic smooth spherocylinders
FC      ?= gfortran
FFLAGS  ?= -O3 -march=native -funroll-loops
DEBUG   = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra

SRC = src/precision.f90 src/npz_io.f90 src/random.f90 src/closure_model.f90 \
      src/particles.f90 src/collisions.f90 src/distributions.f90 src/dsmc0d.f90
OBJ = $(patsubst src/%.f90,build/%.o,$(SRC))

bin/dsmc0d: $(OBJ) | bin
	$(FC) $(FFLAGS) -o $@ $(OBJ)

build/%.o: src/%.f90 | build
	$(FC) $(FFLAGS) -Jbuild -c $< -o $@

debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="$(DEBUG)"

build bin:
	mkdir -p $@

clean:
	rm -rf build bin

.PHONY: clean debug

# module dependencies
build/npz_io.o build/random.o: build/precision.o
build/closure_model.o: build/npz_io.o build/random.o
build/particles.o: build/random.o
build/collisions.o: build/closure_model.o build/particles.o
build/distributions.o: build/particles.o
build/dsmc0d.o: build/collisions.o build/distributions.o
