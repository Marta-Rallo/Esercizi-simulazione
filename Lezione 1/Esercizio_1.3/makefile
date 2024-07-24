CC = g++

prog3 : main.o random.o buffon.o
	g++ main.o random.o buffon.o -o prog3
main.o : main.cpp
	g++ -c main.cpp -o main.o
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o
buffon.o : EsperimentoBuffon.cpp EsperimentoBuffon.h
	g++ -c EsperimentoBuffon.cpp -o buffon.o
clean :
	rm *.o