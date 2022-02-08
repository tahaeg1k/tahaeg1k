#include "BlackScholes.h"
#include "random_singleton.h"

#include <cassert>
#include <cfloat>
#include <iostream>


BlackScholes::BlackScholes()= default;

BlackScholes::~BlackScholes()= default;

double BlackScholes::callOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return spotPrice* exp(-dividendYield*yearsToExpiry) * Random::normalCDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
            - strike * exp(-riskFreeInterestRate * yearsToExpiry) * Random::normalCDF(this->d2(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield));
}

double BlackScholes::callOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	if (0.0 == dividendYield)
		this->OptionsPricingModel::callOptionDelta(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield);

	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return exp(-dividendYield*yearsToExpiry) * Random::normalCDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield));
}

double BlackScholes::callOptionGamma(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
    assert(spotPrice >= 0.0);
    assert(strike >= 0.0);
    if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
    assert(volatility >= 0.0);

    //double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
    //double volatility = _volatility / 100.0;
    //double dividendYield = _dividendYield / 100.0;

    return exp(-dividendYield*yearsToExpiry) * Random::normalPDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) / (spotPrice*volatility*std::sqrt(yearsToExpiry));
}

double BlackScholes::callOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return spotPrice * exp(-dividendYield*yearsToExpiry) * std::sqrt(yearsToExpiry) *
            Random::normalPDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield));
}

double BlackScholes::callOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return dividendYield * spotPrice * exp(-dividendYield*yearsToExpiry) * Random::normalCDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
            -riskFreeInterestRate * strike * exp(-riskFreeInterestRate*yearsToExpiry) * Random::normalCDF(this->d2(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
            - spotPrice * exp(-dividendYield*yearsToExpiry)* Random::normalPDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) * volatility/(2*std::sqrt(yearsToExpiry));
}

double BlackScholes::callOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return yearsToExpiry * strike * std::exp(-riskFreeInterestRate*yearsToExpiry) *
            Random::normalCDF(this->d2(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) ;
}

double BlackScholes::putOptionValue(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return -spotPrice * std::exp(-dividendYield*yearsToExpiry) * Random::normalCDF(-this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
     + exp(-riskFreeInterestRate * yearsToExpiry) * strike * Random::normalCDF(-this->d2(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield));
}

double BlackScholes::putOptionDelta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return -std::exp(-dividendYield*yearsToExpiry) * Random::normalCDF(-this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield));
}

double BlackScholes::putOptionGamma(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
    assert(spotPrice >= 0.0);
    assert(strike >= 0.0);
    if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
    assert(volatility >= 0.0);

    //double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
    //double volatility = _volatility / 100.0;
    //double dividendYield = _dividendYield / 100.0;

    return exp(-dividendYield*yearsToExpiry) * Random::normalPDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) / (spotPrice*volatility*std::sqrt(yearsToExpiry));
}

double BlackScholes::putOptionVega(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;
    return spotPrice * std::exp(-dividendYield*yearsToExpiry) * std::sqrt(yearsToExpiry) *
           Random::normalPDF(this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) * sqrt(yearsToExpiry);
}

double BlackScholes::putOptionTheta(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	return riskFreeInterestRate * strike * exp(-riskFreeInterestRate*yearsToExpiry) * Random::normalCDF(-this->d2(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
           -dividendYield * spotPrice * exp(-dividendYield*yearsToExpiry) * Random::normalCDF(-this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield))
           - strike * exp(-riskFreeInterestRate*yearsToExpiry)* Random::normalPDF(-this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) * volatility/(2*std::sqrt(yearsToExpiry));
}

double BlackScholes::putOptionRho(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

    return yearsToExpiry * spotPrice * std::exp(-dividendYield*yearsToExpiry) *
            Random::normalCDF(-this->d1(spotPrice, strike, yearsToExpiry, riskFreeInterestRate, volatility, dividendYield)) * exp(-riskFreeInterestRate * yearsToExpiry) * yearsToExpiry;
}

double BlackScholes::d1(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	return (1.0 / volatility * std::sqrt(yearsToExpiry)) * (std::log(spotPrice / strike) + (riskFreeInterestRate - dividendYield + 0.5 * volatility * volatility) * yearsToExpiry);
}

double BlackScholes::d2(double spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	return (1.0 / volatility * std::sqrt(yearsToExpiry)) * (log(spotPrice / strike) + (riskFreeInterestRate - dividendYield - 0.5 * volatility * volatility) * yearsToExpiry);
}