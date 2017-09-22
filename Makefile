
CC = gcc
CFLAG = -g -Wall -I
LIBS = -lm -fopenmp

main: main.o functions.o corrcode.o readinput.o
	$(CC) main.o functions.o corrcode.o readinput.o -o lifetime_corr_omp.exe $(LIBS)

main.o: main.c
	$(CC) -c main.c $(LIBS)

corrcode.o: corrcode.c
	$(CC) -c corrcode.c $(LIBS)

functions.o: functions.c
	$(CC) -c functions.c $(LIBS)

readinput.o: readinput.c
	$(CC) -c readinput.c $(LIBS)
clean:
	rm lifetime_corr_omp.exe main.o functions.o corrcode.o readinput.o

