#include "random.h"
#include <iostream>

using namespace std;

int main() {
  int throws = 10000;
  int blocks = 100;

  Random rnd;

  // Calcolo integrale con distribuzione uniforme(Teo media)
  vector<double> Media = MediaBlocchi_RandMedia(blocks, throws, rnd, 0, 1);
  vector<double> Mediablocchi = mediablocchi(blocks, throws, Media);
  vector<double> Mediablocchi2 = mediablocchi2(blocks, throws, Media);
  vector<double> Err = Varianza(blocks, Mediablocchi, Mediablocchi2);

  // Calcolo integrale con importance sampling
  vector<double> MediaHoM =
      MediaBlocchi_RandHoM(blocks, throws, rnd, 0, 1, 3. / 2);
  vector<double> MediablocchiHoM = mediablocchi(blocks, throws, MediaHoM);
  vector<double> Mediablocchi2HoM = mediablocchi2(blocks, throws, MediaHoM);
  vector<double> ErrHoM = Varianza(blocks, MediablocchiHoM, Mediablocchi2HoM);

  WriteonFile("ris.dat", Mediablocchi, Err, MediablocchiHoM, ErrHoM);

  return 0;
}