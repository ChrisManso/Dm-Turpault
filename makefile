exe:fonctions.o main.o
	gfortran -o exe fonctions.o main.o -g -fbounds-check
fonctions.o:fonctions.f90
	gfortran -c fonctions.f90
main.o:main.f90 fonctions.o
	gfortran -c main.f90
