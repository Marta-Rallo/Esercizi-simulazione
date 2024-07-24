CC = g++
CFLAGS = -Wall -O3 --std=c++11

genetic : genetic.o random.o
	$(CC) random.o genetic.o -g -o genetic

genetic.o : genetic.cpp
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out