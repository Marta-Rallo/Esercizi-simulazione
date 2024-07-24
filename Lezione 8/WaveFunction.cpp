#include "WaveFunction.h"

WaveFunction::WaveFunction(double sigma, double mu) : sigma(sigma), mu(mu) {}

double WaveFunction::PsiTrial(double x) const {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma)) + exp(-pow(x + mu, 2) / (2 * sigma * sigma));
}

double WaveFunction::PsiSquared(double x) const {
    double psi = PsiTrial(x);
    return psi * psi;
}

double WaveFunction::Energy(double x) const {
    double sigma2 = sigma * sigma;
    double exp1 = exp(-pow(x - mu, 2) / (2 * sigma2));
    double exp2 = exp(-pow(x + mu, 2) / (2 * sigma2));
    double kinetic = (exp1/sigma2)*(((x-mu)*(x-mu)/sigma2)-1)+(exp2/sigma2)*(((x+mu)*(x+mu)/sigma2)-1);
    double potential = pow(x, 4) - 2.5 * x * x;
    return (-((hbar * hbar) / (2 * m)) * kinetic) / PsiTrial(x) + potential;
}
