#FC = ifort
FC = mpif90
CC = cc
PROG =  exchange
# SUN
#LIBS = -lpthread -libguide -lmkl -openmp -LFFTW/lib -lfftw3 -lm -L/usr/local/mpi_qlogic/lib64
LIBS = -L/usr/local/mpi_qlogic/lib64 -L/usr/local/intel/PStudio_XE/mkl/lib/intel64 -lpthread -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -LFFTW/lib -lfftw3 -lm -L/usr/local/mpi_qlogic/lib64
FFLAGS = -O3
IFLAGS = -IFFTW/include -I/usr/local/mpi_qlogic/include



OBJ = parameters.o functions.o fftw3.o exchange.o 

all: $(PROG) $(PROGO)

clobber:
	rm -f *.o

clean:clobber
	rm -f exchange

$(PROG): $(OBJ)
	$(FC) $(FFLAGS) -o $(PROG) $(OBJ) $(LIBS) -static-intel

#mpimod.o : mpimod.f90 
#	$(FC) $(FFLAGS) $(IFLAGS) -c mpimod.f90

parameters.o : parameters.f90 
	$(FC) $(FFLAGS) -c parameters.f90

functions.o : functions.f90 
	$(FC) $(FFLAGS) -c functions.f90

fftw3.o : fftw3.f90 
	$(FC) $(FFLAGS) $(IFLAGS) -c fftw3.f90

exchange.o : exchange.f90 
	$(FC) $(FFLAGS) -c exchange.f90

.f.o : 
	$(FC) $(FFLAGS) -c  $<

.c.o : 
	$(CC) -O2  -c  $<
