
CC = gcc
CFLAG = -g -Wall -I -O3 
LIBS = -lm 

main: main.o functions.o corrcode.o readinput.o head.h
	$(CC) main.o functions.o corrcode.o readinput.o -o lifetime_corr.exe $(LIBS) $(CFLAG)

main.o: main.c head.h
	$(CC) -c main.c $(CFLAG)

corrcode.o: corrcode.c head.h
	$(CC) -c corrcode.c $(CFLAG)

functions.o: functions.c head.h
	$(CC) -c functions.c $(CFLAG)

readinput.o: readinput.c head.h
	$(CC) -c readinput.c $(CFLAG)

clean:
	rm lifetime_corr.exe main.o functions.o corrcode.o readinput.o

