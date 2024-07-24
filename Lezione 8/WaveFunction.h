#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include "random.h"

class WaveFunction {
public:
    WaveFunction(double sigma, double mu);
    double PsiTrial(double x) const;
    double PsiSquared(double x) const;
    double Energy(double x) const;

private:
    double sigma;
    double mu;
    static constexpr double hbar = 1.0;
    static constexpr double m = 1.0;
};

#endif // WAVEFUNCTION_H
