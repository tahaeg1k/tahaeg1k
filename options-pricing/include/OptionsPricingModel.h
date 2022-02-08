/**
 * @file OptionsPricingModel.h
 * @brief Base class for any of the models or methods used for pricing:
 *      - Black Scholes
 *      - Monte Carlo by Malliavin
 *
 * All functions are virtual and for callOptionValue and putOptionValue they must be
 * implemented in the derived classes.
 *
 * @author BENNANI Mohamed Taha & LAKHDAR Othmane
 *
 * @bugs No known bugs for now.
 */

#pragma once


class OptionsPricingModel
{
public:
	OptionsPricingModel();
	virtual ~OptionsPricingModel();

	virtual double callOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const=0;
	virtual double callOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;

	virtual double putOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const=0;
	virtual double putOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
};
