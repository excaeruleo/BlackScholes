#ifndef EXCBS_SOLVER_H_
#define EXCBS_SOLVER_H_

#include <vector>

namespace excbs {

// Black-Scholes Euler finite-difference solver
class BlackScholesSolver {
public:

    // Analytical Black-Scholes solution for comparison
    static double normalCDF(double x) noexcept;

    // ctor. requires risk-free rate and diffusion constant
    BlackScholesSolver(double risk_free_rate, double volatility) noexcept; 

    // Black-Scholes analytical call price
    double blackScholesCall(double S, double K, double T, double t = 0) const;

    // Numerical solution using finite differences
    std::vector<double> solveNumerical(double S0, double K, double T, 
                                     int numS, int numT) const;

private:
    double r;      // Risk-free rate
    double sigma;  // Volatility
};

}  // namespace excbs

#endif  // EXCBS_SOLVER_H_
