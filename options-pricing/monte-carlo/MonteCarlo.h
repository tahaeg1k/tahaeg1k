/** @file MonteCarlo.h
 *  @brief Monte Carlo simulations using Malliavin calculus.
 *
 *  This file contains the function declarations for Monte Carlo simulations.
 *
 *  @author Mohamed Taha Bennani
 *  @bug No known bugs.
 *
 */

#pragma once

#include "OptionsPricingModel.h"
#include "random_singleton.h"

#include <string>
#include <vector>
#include <map>

/**
 * @brief Class MonteCarlo that contains two function for computing the price
 *        and the Greeks using the Malliavin method.
 *
 * Two types of products :
 *      - European with payoff g(ST)
 *      - Exotic with payoff g( (St){t in [0,T]} ), example Asian call option.
 */

class MonteCarlo : public OptionsPricingModel
{
public:
    const int number_of_iterations= 100000;
    const static int length_path = 100;// Subdivision of T


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
