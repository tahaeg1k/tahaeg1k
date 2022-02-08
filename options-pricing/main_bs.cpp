//
// Created by Mohamed Taha Bennani on 21/01/2022.
//

#include "BlackScholes.h"
#include "MonteCarlo.h"

#include <iostream>

int main(){
    std::cout << "Hello, World!" << std::endl;

    double S = 100.0;  // spot price
    double K = 50.0;  // Strike price
    double T = 1.0;    // One year until expiry
    double r = 0.05;   // Risk-free rate (5%)
    double v = 0.30;    // Volatility of the underlying (30%)
    double q = 0.0; // dividend yield of the stock

    BlackScholes model_bs;
    double call_t_bs = model_bs.callOptionValue(S, K, T, r, v, q);
    std::cout << "Call option value using Black-Scholes formula: " << call_t_bs << std::endl;
    double delta_bs = model_bs.callOptionDelta(S, K, T, r, v, q);
    std::cout << "Delta call  option value using Black-Scholes formula: " << delta_bs << std::endl;

    MonteCarlo model_mc;
    double call_t_mc = model_mc.callOptionValue(S,K,T,r,v,q);
    std::cout << "Call option value using Monte Carlo : " << call_t_mc << std::endl;
    double delta_mc = model_mc.callOptionDelta(S, K, T, r, v, q);
    std::cout << "Delta call option value using finite difference : " << delta_mc << std::endl;

    return 0;

}

