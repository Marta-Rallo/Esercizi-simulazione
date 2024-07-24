CC = g++
CFLAGS = -Wall -O3 --std=c++11

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

es8.1 : es8.1.o random.o
	$(CC) random.o es8.1.o -o es8.1

es8.1.o : es8.1.cpp
	$(CC) -c es8.1.cpp -o es8.1.o $(CFLAGS)

es8.2 : es8.2.o random.o
	$(CC) random.o es8.2.o -o es8.2

es8.2.o : es8.2.cpp
	$(CC) -c es8.2.cpp -o es8.2.o $(CFLAGS)

main : main.o random.o WaveFunction.o SimulatedAnnealing.o
	$(CC) random.o WaveFunction.o SimulatedAnnealing.o main.o -o main

clean :
	rm *.o main.exe seed.out
