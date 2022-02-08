/** @file BlackScholes.h
 *  @brief Function declarations for the Black-Scholes model.
 *
 *  This file contains the function declarations for the pricing and Greeking of the European call and put.
 *  The functions are based on the explicit formulas in the Black-Scholes, for the call and the put.
 *
 *  @author Mohamed Taha Bennani
 *  @bug No known bugs
 */

#pragma once

#include "OptionsPricingModel.h"

/** @class BlackScholes
 *  @brief BlackScholes class, subclass of OptionsPricingModel
 *
 *  For now, computation is based on the cumulative distribution function N
 *  defined in random_singleton.h
 */
class BlackScholes : public OptionsPricingModel
{
public:
	BlackScholes();
	virtual ~BlackScholes();

public:
    /**
     * @brief Call option value
     * @param spotPrice
     * @param strike
     * @param yearsToExpiry
     * @param riskFreeInterestRate
     * @param volatility
     * @param dividendYield
     * @return Call option value as a double.
     */
	virtual double callOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionGamma(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
    virtual double callOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double callOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;

    /**
     * @brief Put option value in the model of Black-Scholes.
     * @param spotPrice
     * @param strike
     * @param yearsToExpiry
     * @param riskFreeInterestRate
     * @param volatility
     * @param dividendYield
     * @return Value of a put option as a double.
     */
	virtual double putOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
    virtual double putOptionGamma(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
    virtual double putOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	virtual double putOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;

private:
	inline double d1(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
	inline double d2(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const;
};
