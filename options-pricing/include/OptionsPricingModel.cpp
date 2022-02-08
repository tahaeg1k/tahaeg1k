#include "OptionsPricingModel.h"

#include <limits>
#include <cstdlib>
#include <iostream>
#include <cfloat>


OptionsPricingModel::OptionsPricingModel() = default;

OptionsPricingModel::~OptionsPricingModel() = default;

double OptionsPricingModel::callOptionDelta(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->callOptionValue(spotPrice + FLT_EPSILON, strike, yearsToExpiry, _riskFreeInterestRate, _volatility, _dividendYield) -
                   this->callOptionValue(spotPrice - FLT_EPSILON, strike, yearsToExpiry, _riskFreeInterestRate, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::callOptionVega(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->callOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate, _volatility + FLT_EPSILON, _dividendYield) -
                   this->callOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate, _volatility - FLT_EPSILON, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::callOptionTheta(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->callOptionValue(spotPrice, strike, yearsToExpiry + FLT_EPSILON, _riskFreeInterestRate, _volatility, _dividendYield) -
                   this->callOptionValue(spotPrice, strike, yearsToExpiry - FLT_EPSILON, _riskFreeInterestRate, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::callOptionRho(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->callOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate + FLT_EPSILON, _volatility, _dividendYield) -
                   this->callOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate - FLT_EPSILON, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::putOptionDelta(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->putOptionValue(spotPrice + FLT_EPSILON, strike, yearsToExpiry, _riskFreeInterestRate, _volatility, _dividendYield) -
                   this->putOptionValue(spotPrice - FLT_EPSILON, strike, yearsToExpiry, _riskFreeInterestRate, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::putOptionVega(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->putOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate, _volatility + FLT_EPSILON, _dividendYield) -
                   this->putOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate, _volatility - FLT_EPSILON, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::putOptionTheta(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->putOptionValue(spotPrice, strike, yearsToExpiry + FLT_EPSILON, _riskFreeInterestRate, _volatility, _dividendYield) -
                   this->putOptionValue(spotPrice, strike, yearsToExpiry - FLT_EPSILON, _riskFreeInterestRate, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}

double OptionsPricingModel::putOptionRho(double spotPrice, double strike, double yearsToExpiry, double _riskFreeInterestRate, double _volatility, double _dividendYield) const {
	return (
                   this->putOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate + FLT_EPSILON, _volatility, _dividendYield) -
                   this->putOptionValue(spotPrice, strike, yearsToExpiry, _riskFreeInterestRate - FLT_EPSILON, _volatility, _dividendYield)
		) / (2.0 * FLT_EPSILON);
}
