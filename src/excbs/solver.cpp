#include "excbs/solver.h"

#include <cmath>

namespace excbs {

BlackScholesSolver::
BlackScholesSolver(double risk_free_rate, double volatility) noexcept 
  : r(risk_free_rate), sigma(volatility)
{}

double 
BlackScholesSolver::normalCDF(double x) noexcept
{
    return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
}

double
BlackScholesSolver::
blackScholesCall(double S, double K, double T, double t) const
{
    double tau = T - t;
    if (tau <= 0) return std::max(S - K, 0.0);

    double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*tau) / (sigma*sqrt(tau));
    double d2 = d1 - sigma*sqrt(tau);

    return S*normalCDF(d1) - K*exp(-r*tau)*normalCDF(d2);
}

std::vector<double>
BlackScholesSolver::
solveNumerical(double S0, double K, double T, int numS, int numT) const
{
    double Smax = 2 * std::max(S0, K);  // Maximum stock price
    double Smin = 0;                    // Minimum stock price
    double dS = (Smax - Smin) / (numS - 1);
    double dt = T / numT;

    // Initialize grid
    std::vector<double> V_old(numS);
    std::vector<double> V_new(numS);
    std::vector<double> S(numS);

    // Set up stock price array and initial condition (t=0)
    for (int i = 0; i < numS; i++) {
        S[i] = i * dS;
        // Initial condition (payoff at t=0)
        V_old[i] = std::max(S[i] - K, 0.0);
    }

    // Time stepping (forward in time from 0 to T)
    for (int t = 0; t < numT; t++) {
        for (int i = 1; i < numS-1; i++) {
            double delta = (V_old[i+1] - V_old[i]) / dS;  // Forward difference
            double gamma = (V_old[i+1] - 2.0*V_old[i] + V_old[i-1]) / (dS * dS);
            
            V_new[i] = V_old[i] + dt * (
                0.5 * sigma * sigma * S[i] * S[i] * gamma +    
                r * S[i] * delta - 
                r * V_old[i]
            );
        }

        // Boundary conditions
        V_new[0] = 0;  // At S = 0
        V_new[numS-1] = Smax - K * exp(-r * (t*dt));  // At S = Smax

        // Update old values
        V_old = V_new;
    }

    return V_new;
}

}  // namespace excbs
