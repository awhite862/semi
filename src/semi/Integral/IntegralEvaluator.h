#ifndef IntegralEvaluator_H
#define IntegralEvaluator_H

namespace Semi {

//S_uv
double CalculateOverlap(double tau, double rho, double kappa, double rho_alpha, double rho_beta, int *aOrbitalType, int *bOrbitalType);

double findRotation();

} //namespace Semi


#endif //IntegralEvaluator_H
