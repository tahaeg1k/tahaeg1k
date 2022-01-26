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

class MonteCarlo : public OptionsPricingModel
{
private:
	struct Configuration {
		double S;
	};

public:
	MonteCarlo();
	virtual ~MonteCarlo();

public:
	virtual double callOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
};
