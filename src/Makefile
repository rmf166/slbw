# GNU Fortran Compiler
F90    = gfortran
LIBS   =
FFLAGS = -O0 -g -fconvert=big-endian -fno-automatic -fbacktrace \
         -ffpe-trap=invalid,zero,overflow,underflow,denormal \
         -fPIC -fcheck=all -fsanitize=address -fsanitize=undefined \
         -Wall -Wextra -fimplicit-none -Wuninitialized -pedantic \
         -Warray-temporaries

slbw.exe : phys_dat.o main.f90
	$(F90) $(FFLAGS) $(LIBS) -o slbw.exe phys_dat.o main.f90

phys_dat.o : phys_dat.f90
	$(F90) $(FFLAGS) $(LIBS) -c phys_dat.f90

clean:
	rm -f slbw.exe *.mod *.o

