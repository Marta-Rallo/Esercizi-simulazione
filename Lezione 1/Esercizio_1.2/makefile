CC = g++

prog2 : main.o random.o
	g++ main.o random.o -o prog2
main.o : maines2.cpp
	g++ -c maines2.cpp -o main.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o
clean :
	rm *.o