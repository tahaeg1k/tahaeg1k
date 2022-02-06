/** @file MonteCarlo.h
 *  @brief Monte Carlo simulations.
 *
 *  This file contains the function declaration for Monte Carlo simulations.
 *
 *  @author Mohamed Taha Bennani
 *  @bug No known bugs.
 *
 */

#pragma once

#include "OptionsPricingModel.h"
#include "EuropeanVanilla.h"
#include "random_singleton.h"

#include <string>
#include <vector>
#include <map>



class MonteCarlo : public OptionsPricingModel
{
public:
    const int number_of_iterations= 1000000;
    const static int m = 100;// Subdivision of T
    double Price[m];

	struct Configuration {
		double S;
	};


public:
	MonteCarlo();
	virtual ~MonteCarlo();


    std::map<std::string, vector<double>> MalliavinEuropeanVanilla(double (*payoff)(double spotPrice, double strike ),
                                                           const int number_iterations=10000, const int length_path=10,
                                                           double spotPrice=100.0, double strike=50, double yearsToExpiry=1.0,
                                                           double riskFreeInterestRate=0.05, double volatility=0.3,
                                                           double dividendYield=0.0);

    std::map<std::string, vector<double>> MalliavinAsianExotic(double (*payoff)(vector<double> vectorPrice, double strike ),
                                                            const int number_iterations=10000, const int length_path=10,
                                                            double spotPrice=100.0, double strike=50, double yearsToExpiry=1.0,
                                                            double riskFreeInterestRate=0.05, double volatility=0.3,
                                                            double dividendYield=0.0);

    virtual double callOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
};
