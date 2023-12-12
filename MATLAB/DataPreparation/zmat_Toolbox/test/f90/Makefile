F90 :=gfortran

all:
	$(F90) -g -Wall -pedantic -c ../../fortran90/zmatlib.f90
	$(F90) -g -Wall -pedantic testzmat.f90 -o testzmat -L../../lib -lzmat -lz
clean:
	-rm -f testzmat *.o *.mod
