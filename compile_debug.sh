gfortran -g3 -fbounds-check -Wall -Wextra -Warray-temporaries -Wconversion -pedantic -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan matrix_exponential.f90 fftpack5.f90 mainvar.f90 WPevolve_adiabatic_version-1.5.f90 -o WPevolve_adiabatic_version-1.5.x

