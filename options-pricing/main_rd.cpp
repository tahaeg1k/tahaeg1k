//
// Created by Mohamed Taha Bennani on 25/01/2022.
//

#include "random_singleton.h"
#include "BlackScholes.h"
#include "MonteCarlo.h"

#include <iostream>
#include <vector>
#include <map>
#include <fstream>


double CallPayoff(double spotPrice, double strike) {
    return (spotPrice>=strike)? spotPrice-strike:0.0 ;
}

double DigitalPayoff(double spotPrice, double strike) {
    return (spotPrice>=strike)? 1.0:0.0 ;
}

double PutPayoff(double spotPrice, double strike) {
    return (strike>=spotPrice)? strike-spotPrice:0.0 ;
}
/*
std::map<std::string, double> MalliavinEuropeanVanilla(const int number_iterations=10000, const int length_path=10,
                                   double spotPrice=100.0, double strike=50, double yearsToExpiry=1.0,
                                   double riskFreeInterestRate=0.05, double volatility=0.3,
                                   double dividendYield=0.0,
                                   double (*payoff)(double spotPrice, double strike ) = CallPayoff){

    Random::Randomize(time(0));

    double T = yearsToExpiry;
    double h = T/double(length_path);
    double S0 = spotPrice;
    double r = riskFreeInterestRate;
    double q =dividendYield;
    double sigma = volatility;
    double K = strike;

    double normale;
    double price_value = 0;
    double delta_value = 0;
    double gamma_value = 0;
    double vega_value = 0;

    for( int i=0; i<number_iterations; i++){
        std::vector<double> W(length_path+1);
        std::vector<double> Price(length_path+1);
        Price[0] = 100.0;

        for(int j=0; j<=length_path-1; j++){
            normale = Random::Gaussian(0,1);
            W[j+1] = W[j] + std::sqrt(h) * normale;
            Price[j+1] = Price[j] * exp( (r -q - 0.5 * sigma * sigma)*h + sigma* (W[j+1]-W[j]) );
            //mean_squared += (ST >= K) ? (ST - K) * (ST - K) : 0.0;
            //W[j+1] = W[j] + std::sqrt(h)*normale;
        }

        double ST = Price[length_path];
        double WT = W[length_path];
        price_value += CallPayoff(ST, K);
        delta_value += CallPayoff(ST, K)*WT/(S0*sigma*T);
        gamma_value += CallPayoff(ST, K)*(WT*WT - sigma*T*WT - T)/std::pow(S0*sigma*T, 2);
        vega_value += CallPayoff(ST, K)*(WT*WT - sigma*T*WT - T)/(sigma*T);

    }

    price_value *= exp(-r* T); // Discount using the constant interest rate r
    delta_value *= exp(-r* T);
    gamma_value *= exp(-r* T);
    vega_value *= exp(-r* T);

    price_value /= double(number_iterations); // Mean over all simulated samples
    delta_value /= double(number_iterations);
    gamma_value /= double(number_iterations);
    vega_value /= double(number_iterations);

    std::map<std::string, double> greeks {{std::string("Price"), price_value},
                                          {std::string("Delta"), delta_value},
                                          {std::string("Gamma"), gamma_value},
                                          {std::string("Vega"), vega_value}};
    return greeks;

}

*/

void vector_to_csv(const vector<double> simulation, string filename) {
    std::ofstream outfile;
    //std::ios_base::app
    outfile.open("Filename_" + filename + ".txt");
    for (int i = 0; i < simulation.size(); i++) {
        outfile << simulation[i] << std::endl;
    }
    outfile.close();
}

int main(){
    std::cout << "Hello, World!" << std::endl;

    // First initialize the seed :
    Random::Randomize(time(0));

    // Set the number of iterations and the length of the path
    int number_iterations = 10000; // Number of iterations
    const int length_path = 100;
    double T = 1.0;
    double h = T/double(length_path);

    double normale;

    double S0 = 100.0;
    double r = 0.05;
    double q =0.0;
    double sigma = 0.3;
    double K = 50;

    MonteCarlo MC_Call;
    std::map<std::string,vector<double>> greeks_map = MC_Call.MalliavinEuropeanVanilla(CallPayoff);

    double price_value = greeks_map["Price"][number_iterations-1];
    double delta_value = greeks_map["Delta"][number_iterations-1];
    double gamma_value = greeks_map["Gamma"][number_iterations-1];
    double vega_value = greeks_map["Vega"][number_iterations-1];

    // Export vector to Python
    vector_to_csv(greeks_map["Price"], "Price");
    vector_to_csv(greeks_map["Delta"], "Delta");
    vector_to_csv(greeks_map["Gamma"], "Gamma");
    vector_to_csv(greeks_map["Vega"], "Vega");


    // Black Scholes model values
    BlackScholes model_bs;
    double call_price_bs = model_bs.callOptionValue(S0, K, T, r, sigma, 0);
    double call_delta_bs = model_bs.callOptionDelta(S0, K, T, r, sigma, 0);
    double call_vega_bs = model_bs.callOptionVega(S0, K, T, r, sigma, 0);
    double call_gamma_bs = model_bs.callOptionGamma(S0, K, T, r, sigma, 0);



    cout << "------PRICE COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << price_value <<  endl;
    cout << " Theoretical price from Random file= " << Random::Call_price_BS(S0, r, T, K, sigma) << endl;
    std::cout << "Call option value using Black-Scholes formula: " << call_price_bs << std::endl;
    cout << endl;

    cout << "------DELTA COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo Malliavin = " << delta_value <<endl;
    cout << "Call delta value using Black-Scholes formula = " << call_delta_bs <<endl;
    cout << endl;

    cout << "------GAMMA COMPARISON------"<< endl;
    cout << "Empirical gamma found by Monte Carlo = " << gamma_value <<endl;
    cout << "Call gamma value using Black-Scholes formula = " << call_gamma_bs <<endl;
    cout << endl;

    cout << "------VEGA COMPARISON------"<< endl;
    cout << "Empirical vega found by Monte Carlo = " << vega_value <<endl;
    cout << "Call vega value using Black-Scholes formula =  "<< call_vega_bs <<endl;

    return 0;

}
