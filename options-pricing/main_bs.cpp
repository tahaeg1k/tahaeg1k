//
// Created by Mohamed Taha Bennani on 21/01/2022.
//

#include "BlackScholes.h"
#include "Binomial.h"
#include "MonteCarlo.h"

#include <iostream>

int main(){
    std::cout << "Hello, World!" << std::endl;

    double S = 100.0;  // spot price
    double K = 50.0;  // Strike price
    double T = 1.0;    // One year until expiry
    double r = 5.0;   // Risk-free rate (5%)
    double v = 20.0;    // Volatility of the underlying (20%)
    double q= 0.0; // dividend yield of the stock

    BlackScholes model_bs;
    double call_t_bs = model_bs.callOptionValue(S, K, T, r, v, q);
    std::cout << "Call option value using Black-Scholes formula: " << call_t_bs << std::endl;

    Binomial model_b;
    double call_t_b = model_b.callOptionValue(S,K,T,r,v,q);
    std::cout << "Call option value using binomial method : " << call_t_b << std::endl;

    MonteCarlo model_mc;
    double call_t_mc = model_mc.callOptionValue(S,K,T,r,v,q);
    std::cout << "Call option value using Monte Carlo : " << call_t_mc << std::endl;

    return 0;

}

