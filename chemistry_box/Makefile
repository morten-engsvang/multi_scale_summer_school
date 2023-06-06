OPT = -O0 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none

all: chemistry_box.exe

chemistry_box.exe: chemistry_box.f90 chemistry_mod.o opkdmain.o opkda1.o opkda2.o
	gfortran $(OPT) chemistry_box.f90 chemistry_mod.o opkdmain.o opkda1.o opkda2.o -o chemistry_box.exe

chemistry_mod.o: chemistry_mod.f90 opkdmain.o opkda1.o opkda2.o
	gfortran $(OPT) -c chemistry_mod.f90

opkda1.o: opkda1.f
	gfortran $(OPT) -std=legacy -w -c opkda1.f

opkda2.o: opkda2.f
	gfortran $(OPT) -std=legacy -c opkda2.f

opkdmain.o: opkdmain.f
	gfortran $(OPT) -std=legacy -c opkdmain.f

clean:
	rm -f *.exe
	rm -f *.o
	rm -f *.mod
