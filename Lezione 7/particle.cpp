/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <math.h>
#include "particle.h"

using namespace std;

void Particle :: initialize(){ //Metto spin a 1 e inizializzo i vettori a una dimensione specifica (ndim)
   _spin = 1;
   _x.resize(_ndim);
   _xold.resize(_ndim);
   _v.resize(_ndim);
   return;
}

void Particle :: translate(vec delta, vec side){ 

//Sposta la posizione di una particella, delta è la quantità di spostamento per ciascuna dimensione, size la dimensione del sistema (scatola),
//e usa le pbc per mantenere la particella nel sistema

   for(unsigned int i=0; i<_ndim; i++){
     _x(i) = pbc(_x(i) + delta(i), side(i));
   }
}

void Particle :: flip(){ //Inverto lo spin della particella
   _spin = -1*this->getspin();
}

void Particle :: moveback(){ 

//ripristina la posizione della particella alla sua posizione precedente (memorizzata in _xold). 
//È utilizzata quando un'operazione di movimento della particella non è accettabile e si desidera tornare allo stato precedente
   _x = _xold;
}

void Particle :: acceptmove(){

//Aggiorna la posizione precedente della particella con la sua posizione attuale. 
//È utilizzata quando un'operazione di movimento della particella è accettabile e si desidera memorizzare la nuova posizione come posizione precedente.
   _xold = _x;
}

int Particle :: getspin(){ //restituisce lo spin della particella
   return _spin;
}

void Particle :: setspin(int spin){ //imposta lo spin della particella col valore specificato
   _spin = spin;
   return;
}

double Particle :: getposition(int dim, bool xnew){ //Serve per ottenere la posizione di una particella lungo una dimensione specificata
   if(xnew) return _x(dim); //Restituisce la posizione corrente lungo la dimensione specificata
   else return _xold(dim); //Restituisce la posizione precedente lungo la dimensione specificata
}

void Particle :: setposition(int dim, double position){ //Imposta la posizione della particella lungo la dimensione specificata
   _x(dim) = position;
   return;
}

void Particle :: setpositold(int dim, double position){ //Imposta la posizione precedente della particella lungo la dimensione specificata
   _xold(dim) = position;
   return;
}

double Particle :: getvelocity(int dim){ //Restituisce la velocità della particella lungo la dimensione specificata
   return _v(dim);
}

vec Particle :: getvelocity(){ //Restituisce un vettore che rappresenta la velocità della particella in tutte le dimensioni
   return _v;
}

void Particle :: setvelocity(int dim, double velocity){ //Imposta la velocità della particella lungo la dimensione specificata
   _v(dim) = velocity;
   return;
}

double Particle :: pbc(double position, double side){ //utilizzata per applicare le pdb ad una posizione specifica, garantendo che la posizione rimanga all'interno dei limiti del sistema simulato
  
  return position - side * rint(position / side); 
  
  //rint() è una funzione che restituisce il valore intero più vicino 
  //(utilizzata per garantire che la posizione finale dopo l'applicazione delle PBC sia un multiplo intero della dimensione del sistema)
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
