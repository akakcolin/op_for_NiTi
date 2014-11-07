CC=gcc
CCFLAGS = -Wall -O2 -g -lm
OMPFLAGS = -fopenmp
.PHONY : clean

all: serial parallel

serial: serial.o 
	$(CC) $(CCFLAGS) $^ -o $@

serial.o: op_serial.c
	$(CC) $(CCFLAGS) $^ -c -o $@

parallel: parallel.o
	$(CC) $(OMPFLAGS) $(CCFLAGS) $^ -o $@

parallel.o: op_parallel.c
	$(CC) $(OMPFLAGS) $(CCFLAGS) $^ -c -o $@

clean:
	rm -rf *.o
