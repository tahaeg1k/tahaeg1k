//
// Created by Mohamed Taha Bennani on 25/01/2022.
//

#include "random_singleton.h"
#include "BlackScholes.h"

#include <iostream>

int main(){
    std::cout << "Hello, World!" << std::endl;

    Random::Randomize(time(0));

    int ite = 1000000; // Number of iterations
    double T = 1;
    double h = T/double(ite);

    double price_value = 0;
    double delta_value = 0;
    double gamma_value = 0;
    double vega_value = 0;
    double var_quad = 0;

    double normale;
    double * moy_courante = new double[ite];
    double number_realis = 0;

    int ite_ext = 10000;
    double ST;
    double S0 = 100;
    double r = 0.05;
    double sigma = 0.3;
    double K = 50;
    double mean_squared = 0;
    double WT;
    double W[ite];
    W[0] = 0;

    for( int j=0; j<ite; j++)
    {
        normale = Random::Gaussian(0,1);
        WT = std::sqrt(T) * normale;
        ST = S0 * exp( (r - 0.5 * sigma * sigma)*T + sigma* WT );
        price_value += (ST >= K) ? ST - K : 0.0;
        delta_value += (ST >= K) ? (ST - K)*WT/(S0*sigma*T) : 0.0;
        gamma_value += (ST >= K) ? (ST - K)*(WT*WT - sigma*T*WT - T)/std::pow(S0*sigma*T, 2) : 0.0;
        vega_value += (ST >= K) ? (ST - K)*(WT*WT - sigma*T*WT - T)/(sigma*T) : 0.0;
        mean_squared += (ST >= K) ? (ST - K) * (ST - K) : 0.0;
        W[j+1] = W[j] + std::sqrt(h)*normale;
    }

    price_value *= exp(-r* T); // Discount using the constant interest rate r
    delta_value *= exp(-r* T);
    gamma_value *= exp(-r* T);
    vega_value *= exp(-r* T);

    price_value /= double(ite); // Mean over all simulated samples
    delta_value /= double(ite);
    gamma_value /= double(ite);
    vega_value /= double(ite);

    BlackScholes model_bs;
    double call_price_bs = model_bs.callOptionValue(S0, K, T, r, sigma, 0);
    double call_delta_bs = model_bs.callOptionDelta(S0, K, T, r, sigma, 0);
    double call_vega_bs = model_bs.callOptionVega(S0, K, T, r, sigma, 0);





    cout << "------PRICE COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << price_value <<  endl;
    cout << " Theoretical price from Random file= " << Random::Call_price_BS(S0, r, T, K, sigma) << endl;
    std::cout << "Call option value using Black-Scholes formula: " << call_price_bs << std::endl;
    cout << endl;

    cout << "------DELTA COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << delta_value <<endl;
    cout << "Call option value using Black-Scholes formula: " << call_delta_bs <<endl;
    cout << endl;

    cout << "------GAMMA COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << gamma_value <<endl;
    cout << "Call option value using Black-Scholes formula: " << "NO FUNCTION IMPLEMENTED YET" <<endl;
    cout << endl;

    cout << "------VEGA COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << vega_value <<endl;
    cout << "Call option value using Black-Scholes formula: "<< call_vega_bs <<endl;

    return 0;

}
