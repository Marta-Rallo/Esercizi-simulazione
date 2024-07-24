CC = g++

	
prog1 : main.o random.o
	g++ main.o random.o -o prog1
	
main.o : main.cpp
	g++ -c main.cpp -o main.o
	
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o
	
clean :
	rm *.o