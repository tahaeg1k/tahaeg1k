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
    int Bin;
    int N = 10;
    double Ber;
    double prob = 0.5;
    double expo;
    double cauchy;
    double lambda = 1.;
    double gamma;
    double a = 0.5;
    double T = 1;
    double h = T/double(ite);

    double mean = 0;
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

    double W[ite];

    for( int j=0; j<ite-1; j++)
    {
        //  Bin = Random::Binomial(prob);
        // Ber = Random::Bernoulli(N,prob);
        // mean += Bin;
        // expo = Random::Exponentielle(lambda);
        // cauchy = Random::Cauchy(0, lambda);
        // normale = Random::BoxMuller(0, 1);
        normale = Random::Gaussian(0,1);
        ST = S0 * exp( (r - 0.5 * sigma * sigma)*T + sigma * sqrt(T) * normale);
        mean += (ST >= K) ? ST - K : 0.0;
        mean_squared += (ST >= K) ? (ST - K) * (ST - K) : 0.0;
        // var_quad += normale * normale;
        //  mean += expo;
        // moy_courante[j+1] = mean/(double(j+1));
        W[j+1] = W[j] + normale;
        // gamma = Random::Gamma(a);
        // mean += gamma;
        // cout << moy_courante[j+1] << endl;
        // mean  = W[j+1]/((j+1)*h);
    }

    mean /= double(ite); // Mean over all simulated samples
    mean *= exp(-r* T); // Discount using the constant interest rate r
    cout << "Empirical mean found by Monte Carlo = " << mean <<  endl;
    cout << " Theoretical price from Random file= " << Random::Call_price_BS(S0, r, T, K, sigma) << endl;

    BlackScholes model_bs;
    double call_t_bs = model_bs.callOptionValue(S0, K, T, r, sigma, 0);
    std::cout << "Call option value using Black-Scholes formula: " << call_t_bs << std::endl;


    return 0;

}
