#define _USE_MATH_DEFINES


#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>

// These libraries provide access to cout, math functions such as exponential and log, and max function.
// First function to be implemented is the Box-Muller algorithm, which allows to generate a Gaussian
// random variable from two uniform random variables. Later on, we will use std::normal_distribution<>
// template class found in the <random> library in order to avoid reinventing the wheel.

double gaussian_box_muller(){
    double x = 0.0;
    double y = 0.0;
    double euclid_sq = 0.0;

    /* initialize random seed: */
    //srand (time(NULL));

    // Continue generating two uniform random variables until the square of their euclidean distance is
    // less than one. This is the polar Box-Muller method.
    do{
        x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
        y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
        euclid_sq = x*x + y*y;
    }while(euclid_sq >= 1.0);
    // This is not the only possible method.
    return x*sqrt( -2*log(euclid_sq) / euclid_sq);
}

// Standard normal probability density function
double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// An approximation to the cumulative distribution function
// for the standard normal distribution
// Note: This is a recursive function.
// More on simulating the norm_cdf function later.
double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

// This calculates d_j, for j in {1,2}. This term appears in the closed
// form solution for the European call or put price
double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}

// Calculate the European vanilla call price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return S * norm_cdf(d_j(1, S, K, r, v, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, v, T));
}

// Calculate the European vanilla put price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return -S*norm_cdf(-d_j(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, v, T));
}


// Pricing a European vanilla call option with a Monte Carlo method
double monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r,
                              const double& v, const double& T){
    /*
     * This function is for pricing a European vanilla call option with a Monte Carlo method.
     *
     * @param num_sins is the number of simulation
     * @param S is for the spot price
     * @param K is the strike
     * @param r is for the interest rate
     * @param v is for the volatility
     * @param T is the maturity
     * @return call option price using the Monte Carlo method.
     */

    double S_adjust = S * exp(T*(r-0.5*v*v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i=0; i<num_sims; i++){
        double gauss_bm = gaussian_box_muller(); // draw of the gaussian law
        S_cur = S_adjust * exp(v*sqrt(T) * gauss_bm);
        payoff_sum += std::max(S_cur - K, 0.0);
    }
    return exp(-r*T) * (payoff_sum / static_cast<double>(num_sims)) ;
}

// Pricing a European vanilla call option with a Monte Carlo method
double monte_carlo_put_price(const int& num_sims, const double& S, const double& K, const double& r,
                              const double& v, const double& T){
    /*
     * This function is for pricing a European vanilla put option with a Monte Carlo method.
     *
     * @param num_sins is the number of simulation
     * @param S is for the spot price
     * @param K is the strike
     * @param r is for the interest rate
     * @param v is for the volatility
     * @param T is the maturity
     * @return call option price using the Monte Carlo method.
     */

    double S_adjust = S * exp(T*(r-0.5*v*v));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    for (int i=0; i<num_sims; i++){
        double gauss_bm = gaussian_box_muller(); // draw of the gaussian law
        S_cur = S_adjust * exp(v*sqrt(T) * gauss_bm);
        payoff_sum += std::max(K - S_cur, 0.0);
    }
    return exp(-r*T) * (payoff_sum / static_cast<double>(num_sims)) ;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    // Let's test the Gaussian box muller method.
    double gauss_bm = gaussian_box_muller();
    std::cout << "This is a draw from the normal distribution law : " << gauss_bm << std::endl;

    int num_sims = 10000000;   // Number of simulated asset paths
    double S = 100.0;  // Option price
    double K = 100.0;  // Strike price
    double r = 0.05;   // Risk-free rate (5%)
    double v = 0.2;    // Volatility of the underlying (20%)
    double T = 1.0;    // One year until expiry
    double mt_call_price = monte_carlo_call_price(num_sims, S ,K, r, v, T);
    std::cout << "This is the result of the Monte Carlo simulation for the call price "
                 ": " << mt_call_price << std::endl;
    double mt_put_price = monte_carlo_put_price(num_sims, S ,K, r, v, T);
    std::cout << "This is the result of the Monte Carlo simulation for the put price "
                 ": " << mt_put_price << std::endl;
    // Then we calculate the call/put analytical values
    double call_analytical = call_price(S, K, r, v, T);
    double put_analytical = put_price(S, K, r, v, T);

    std::cout << "Analytical Call Price:      " << call_analytical << std::endl;
    std::cout << "Analytical Put Price:       " << put_analytical << std::endl;

    return 0;
}
