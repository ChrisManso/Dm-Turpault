#Objet
OBJ= fonctions.o fonctionCSR.o main.o
optn= -g -fbounds-check -o0 -Wall


#Compile et crée l'executable
exe: $(OBJ)
	gfortran $(optn) -o exe $(OBJ)

%.o: %.f90
	gfortran -c $<

#Cleaner
clean:
	rm -r *~ *.o *.mod
