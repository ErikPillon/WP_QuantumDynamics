FC=gfortran
FFLAGS=-w -fallow-argument-mismatch
SRC=jacobi.f90 matrix_exponential.f90 fftpack5.f90 mainvar.f90 mainmethods.f90 WPevolve_adiabatic_version-1.5.f90
OBJ=${SRC:.f90=.o} 

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

WPevolve.x: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: 
	@rm -f *.mod *.o 
