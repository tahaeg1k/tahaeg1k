/**
 * @file main_asian.cpp
 * @brief Main file for Monte Carlo simulations applied to Asian options.
 */

#include "random_singleton.h"
#include "MonteCarlo.h"

#include <iostream>
#include <numeric>
#include <vector>
#include <map>
#include <fstream>


void vector_to_csv(const vector<double> simulation, string filename) {
    std::ofstream outfile;
    //std::ios_base::app
    outfile.open("File_" + filename + ".txt");
    for (int i = 0; i < simulation.size(); i++) {
        outfile << simulation[i] << std::endl;
    }
    outfile.close();
}


double AsianCallPayoff(vector<double> vectorPrice, double strike) {
    double average = std::accumulate(vectorPrice.begin(), vectorPrice.end(), 0.0)
                                        / vectorPrice.size();
    return (average>=strike)? average-strike:0.0 ;
}


int main(){
    std::cout << "Hello, World!" << std::endl;

    // First initialize the seed :
    Random::Randomize(time(0));

    // Set the number of iterations and the length of the path
    int number_iterations = 10000; // Number of iterations
    const int length_path = 1000;
    double T = 1.0;
    double h = T/double(length_path);

    double normale;

    double S0 = 100.0;
    double r = 0.05;
    double q =0.0;
    double sigma = 0.3;
    double K = 50;

    MonteCarlo MC_AsianCall;
    std::map<std::string,vector<double>> greeks_map = MC_AsianCall.MalliavinAsianExotic(AsianCallPayoff,
                                                                                        number_iterations, length_path);

    double price_value = greeks_map["Price"][number_iterations-1];
    double delta_value = greeks_map["Delta"][number_iterations-1];

    // Export vector to Python
    vector_to_csv(greeks_map["Price"], "Asian_Price");
    vector_to_csv(greeks_map["Delta"], "Asian_Delta");


    cout << "------PRICE COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo = " << price_value <<  endl;
    cout << endl;

    cout << "------DELTA COMPARISON------"<< endl;
    cout << "Empirical mean found by Monte Carlo Malliavin = " << delta_value <<endl;
    cout << endl;

    return 0;

}
