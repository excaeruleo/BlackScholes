#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

class BlackScholesSolver {
private:
    double r;      // Risk-free rate
    double sigma;  // Volatility
    
public:
    BlackScholesSolver(double risk_free_rate, double volatility) 
        : r(risk_free_rate), sigma(volatility) {}

    // Analytical Black-Scholes solution for comparison
    double normalCDF(double x) {
        return 0.5 * (1 + erf(x / sqrt(2.0)));
    }

    double blackScholesCall(double S, double K, double T, double t = 0) {
        double tau = T - t;
        if (tau <= 0) return std::max(S - K, 0.0);

        double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*tau) / (sigma*sqrt(tau));
        double d2 = d1 - sigma*sqrt(tau);

        return S*normalCDF(d1) - K*exp(-r*tau)*normalCDF(d2);
    }

    // Numerical solution using finite differences
    std::vector<double> solveNumerical(double S0, double K, double T, 
                                     int numS, int numT) {
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
};

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

    BlackScholesSolver solver(r, sigma);
    auto numerical = solver.solveNumerical(S0, K, T, numS, numT);
    // Print results and compare with analytical solution
    std::cout << "Stock Price\tNumerical\tAnalytical\tDifference\n";
    double dS = (2*std::max(S0, K)) / (numS - 1);
    
    for (int i = 0; i < numS; i += 5) {
        double S = i * dS;
        double analytical = solver.blackScholesCall(S, K, T);
        std::cout << S << "\t\t" 
                 << printf("%.5f", numerical[i]) << "\t"
                 << printf("%.5f", analytical) << "\t"
                 << printf("%.5f", std::abs(numerical[i] - analytical)) << "\n";
    }

    return 0;
} 
