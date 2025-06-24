
all: 
	gfortran -o run.out -fopenmp lapack_mods.f90 algebra.f90 Hamiltonian.f90 mod_Gfunctions.f90 stypesub.f90 StypeJunc.f90 -L/usr/local/opt/openblas/lib -lopenblas 

clear: 

	rm -f sf* *.o *.mod *.out *~ iteration* fort* err* Volt_Current* Print.dat SpinJunction.*
