build: homework
homework: homework.c
	mpicc -o homework homework.c -lm -Wall -g
serial: homework
	mpirun -np 1 homework imagini.in
distrib: homework
	mpirun -np 4 homework imagini.in
clean:
	rm -f homework
