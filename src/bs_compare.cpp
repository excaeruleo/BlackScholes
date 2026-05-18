#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>

#include "excbs/solver.h"

int main() {
    // Parameters
    double S0 = 100.0;    // Initial stock price
    double K = 100.0;     // Strike price
    double T = 1.0;       // Time to maturity
    double r = 0.05;      // Risk-free rate
    double sigma = 0.2;   // Volatility
    
    // Grid parameters
    int numS = 100;       // Number of stock price points
    int numT = 1000;      // Number of time points

    excbs::BlackScholesSolver solver(r, sigma);
    auto numerical = solver.solveNumerical(S0, K, T, numS, numT);
    // Print results and compare with analytical solution
    std::cout <<
        "Stock Price     "
        "Numerical       "
        "Analytical      "
        "Difference\n";
    double dS = (2*std::max(S0, K)) / (numS - 1);
    // column width when printing
    constexpr int col_width = 16;
    // ensure left-aligned output
    std::cout << std::left;
    // line formatting helper for fixed column width + precision
    auto fmt = []() -> auto&
    {
        return std::cout << std::setw(col_width) << std::fixed <<
            std::setprecision(6);
    };
    // print values as table
    for (int i = 0; i < numS; i += 5) {
        double S = i * dS;
        double analytical = solver.blackScholesCall(S, K, T);
        fmt() << S;
        fmt() << numerical[i];
        fmt() << analytical;
        fmt() << (std::abs(numerical[i] - analytical)) << "\n";
    }
    // final flush
    std::cout << std::flush;
    return 0;
} 
