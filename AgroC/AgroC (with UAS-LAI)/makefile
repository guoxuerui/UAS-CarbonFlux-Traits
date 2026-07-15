# makefile for soilco2

PGM = agroC
sources = carbon.f90 input.f90 material.f90 output.f90 sink.f90 soilco2.f90 temper.f90 solute.f90 \
          time.f90 watflow.f90 datatypes.f90 geometry.f90 plants.f90 phosphorus.f90 variables.f90 \
          nitrogen.f90 

objects = ${sources:.f90=.o}

SYSTEM = $(shell uname)
# linux
#COMP = gfortran -Wall -g -fbounds-check
ifeq ($(FC),)
FC=$(F90)
endif
ifeq ($(FC),)
FC=gfortran
endif
ifneq ($(SYSTEM),Linux)
# apple
FC = gfortran
endif
DEBUG = "on"
ifeq ($(DEBUG),"on")
#COMP = $(FC) -Wall -g -fbacktrace -fmax-errors=100 -fcheck=all -ffpe-trap=zero,invalid,overflow -ffpe-summary=all #,inexact,denormal,underflow -ffpe-summary=none/all
COMP = $(FC) -ffpe-trap=zero,overflow
LINK = $(FC)
else
#COMP = $(FC) -Ofast
#COMP = $(FC) -fopenmp -Ofast -ftree-parallelize-loops=8 
COMP = $(FC) -Wall -O2
# -ffast-math speeds up about 2.5%
# -Ofast gives different resuts, no convergence in time loop
#LINK = gfortran -static
LINK = $(FC) -fopenmp
endif

%.o: %.f90
	${COMP} -c $<

$(PGM): $(objects)
	$(LINK) -o $@ $(objects)

clean:
	rm -f $(objects) $(PGM) *.mod

run:
	./$(PGM) > ausgabe.txt

tags:
	etags -l fortran $(sources)

sys:
	echo $(SYSTEM)

depend: 
	gfortran -cpp -MM *.f90

print:
	echo $(objects)

time.o : datatypes.o
sink.o: datatypes.o geometry.o variables.o carbon.o
temper.o: datatypes.o
carbon.o: datatypes.o geometry.o material.o variables.o time.o
output.o: datatypes.o geometry.o carbon.o variables.o material.o nitrogen.o solute.o plants.o phosphorus.o
geometry.o: datatypes.o
plants.o: datatypes.o time.o geometry.o time.o material.o
variables.o : datatypes.o 
input.o: datatypes.o geometry.o material.o variables.o time.o plants.o carbon.o output.o solute.o
soilco2.o: datatypes.o geometry.o material.o carbon.o input.o output.o variables.o time.o watflow.o sink.o temper.o plants.o nitrogen.o solute.o phosphorus.o
watflow.o: datatypes.o geometry.o material.o variables.o time.o
material.o: datatypes.o variables.o
nitrogen.o: datatypes.o geometry.o carbon.o solute.o material.o plants.o
solute.o: datatypes.o material.o geometry.o
phosphorus.o: datatypes.o material.o geometry.o carbon.o
