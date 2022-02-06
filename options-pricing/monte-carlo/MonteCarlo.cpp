#include "MonteCarlo.h"
#include "Cdf.h"

#include <cassert>
#include <cfloat>
#include <numeric>




MonteCarlo::MonteCarlo()
{
}

MonteCarlo::~MonteCarlo()
{
}

std::map<std::string, vector<double>> MonteCarlo::MalliavinEuropeanVanilla(double (*payoff)(double spotPrice, double strike ),
                                                                   const int number_iterations, const int length_path,
                                                       double spotPrice, double strike, double yearsToExpiry,
                                                       double riskFreeInterestRate, double volatility,
                                                       double dividendYield){

    Random::Randomize(time(0));

    double T = yearsToExpiry;
    double h = T/double(length_path);
    double S0 = spotPrice;
    double r = riskFreeInterestRate;
    double q =dividendYield;
    double sigma = volatility;
    double K = strike;

    double normale;

    std::vector<double> Price(number_iterations);
    std::vector<double> Delta(number_iterations);
    std::vector<double> Gamma(number_iterations);
    std::vector<double> Vega(number_iterations);

    double price_value = 0;

    double delta_value = 0;
    double gamma_value = 0;
    double vega_value = 0;

    for( int i=0; i<number_iterations; i++){
        std::vector<double> W(length_path+1);
        std::vector<double> S(length_path+1);
        S[0] = S0;

        for(int j=0; j<=length_path-1; j++){
            normale = Random::Gaussian(0,1);
            W[j+1] = W[j] + std::sqrt(h) * normale;
            S[j+1] = S[j] * exp( (r -q - 0.5 * sigma * sigma)*h + sigma* (W[j+1]-W[j]) );

        }

        double ST = S[length_path];
        double WT = W[length_path];
        //Price
        price_value += (*payoff)(ST, K);
        Price[i] = exp(-r* T)*price_value/(double(i+1));
        //Delta
        delta_value += (*payoff)(ST, K)*WT/(S0*sigma*T);
        Delta[i] = exp(-r* T)*delta_value/(double(i+1));

        //Vega
        vega_value += (*payoff)(ST, K)*(WT*WT - sigma*T*WT - T)/(sigma*T);
        Vega[i] = exp(-r* T)*vega_value/(double(i+1));

        //gamma_value += CallPayoff(ST, K)*(WT*WT - sigma*T*WT - T)/std::pow(S0*sigma*T, 2);
        //gamma_value = vega_value/(std::pow(S0, 2) * sigma * T) ;
        Gamma[i] = Vega[i]/(std::pow(S0, 2) * sigma * T) ;


    }

    price_value *= exp(-r* T); // Discount using the constant interest rate r
    delta_value *= exp(-r* T);
    gamma_value *= exp(-r* T);
    vega_value *= exp(-r* T);

    price_value /= double(number_iterations); // Mean over all simulated samples
    delta_value /= double(number_iterations);
    gamma_value /= double(number_iterations);
    vega_value /= double(number_iterations);

    std::map<std::string, vector<double>> greeks {{std::string("Price"), Price},
                                          {std::string("Delta"), Delta},
                                          {std::string("Gamma"), Gamma},
                                          {std::string("Vega"), Vega}};
    return greeks;

}


std::map<std::string, vector<double>> MonteCarlo::MalliavinAsianExotic(double (*payoff)(vector<double> vectorPrice, double strike ),
                                                                           const int number_iterations, const int length_path,
                                                                           double spotPrice, double strike, double yearsToExpiry,
                                                                           double riskFreeInterestRate, double volatility,
                                                                           double dividendYield){

    Random::Randomize(time(0));

    double T = yearsToExpiry;
    double h = T/double(length_path);
    double S0 = spotPrice;
    double r = riskFreeInterestRate;
    double q =dividendYield;
    double sigma = volatility;
    double K = strike;

    double normale;

    std::vector<double> Price(number_iterations);
    std::vector<double> Delta(number_iterations);

    double price_value = 0;

    double delta_value = 0;

    for( int i=0; i<number_iterations; i++){
        std::vector<double> W(length_path+1);
        std::vector<double> S(length_path+1);
        double integrale_stochastique_S = 0.0;
        S[0] = S0;

        for(int j=0; j<=length_path-1; j++){
            normale = Random::Gaussian(0,1);
            W[j+1] = W[j] + std::sqrt(h) * normale;
            S[j+1] = S[j] * exp( (r -q - 0.5 * sigma * sigma)*h + sigma* (W[j+1]-W[j]) );
            integrale_stochastique_S +=  S[j+1] * (W[j+1] - W[j]);
        }

        double mean_S = std::accumulate(S.begin(), S.end(), 0.0)/ S.size();

        double ST = S[length_path];
        double WT = W[length_path];
        double value_payoff = (*payoff)(S, K);
        //Price
        price_value += value_payoff;
        Price[i] = exp(-r* T)*price_value/(double(i+1));
        //Delta
        delta_value += value_payoff*( (2*integrale_stochastique_S)/(S0*sigma*mean_S) + (1/S0));
        Delta[i] = exp(-r* T)*delta_value/(double(i+1));

    }

    price_value *= exp(-r* T); // Discount using the constant interest rate r
    delta_value *= exp(-r* T);

    price_value /= double(number_iterations); // Mean over all simulated samples
    delta_value /= double(number_iterations);

    std::map<std::string, vector<double>> greeks {{std::string("Price"), Price},
                                                  {std::string("Delta"), Delta} };
    return greeks;

}


double MonteCarlo::callOptionValue(double _spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(_spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	const int m = 100000;
	const int n = 100;

	double deltaT = yearsToExpiry / (n - 1); // yearsToExpiry is (T-t)
	double sqrtDeltaT = sqrt(deltaT);
	double spotPrice = _spotPrice * exp((riskFreeInterestRate - dividendYield) * yearsToExpiry);

	Configuration configuration[n];
	double ret = 0.0;
	for (int i = 0; i < m; ++i) {
		configuration[0].S = spotPrice;
		for (int j = 1; j < n; ++j) {
			double X1 = (double)(rand() % 2 * 2 - 1);
			configuration[j].S = configuration[j - 1].S * exp(volatility * X1 * sqrtDeltaT);
		}

		double V = configuration[n - 1].S - strike;
		if (V < 0.0) V = 0.0;

		ret += V;
	}
	ret /= (double)m;
	ret *= exp(-riskFreeInterestRate * yearsToExpiry);

	return ret;
}

double MonteCarlo::putOptionValue(double _spotPrice, double strike, double yearsToExpiry, double riskFreeInterestRate, double volatility, double dividendYield) const {
	assert(_spotPrice >= 0.0);
	assert(strike >= 0.0);
	if (yearsToExpiry < 0.0) yearsToExpiry = 0.0;
	assert(volatility >= 0.0);

	//double riskFreeInterestRate = _riskFreeInterestRate / 100.0;
	//double volatility = _volatility / 100.0;
	//double dividendYield = _dividendYield / 100.0;

	const int m = 100000;
	const int n = 100;

	double deltaT = yearsToExpiry / (n - 1);
	double sqrtDeltaT = sqrt(deltaT);
	double spotPrice = _spotPrice * exp((riskFreeInterestRate - dividendYield) * yearsToExpiry);

	Configuration configuration[n];
	double ret = 0.0;
	for (int i = 0; i < m; ++i) {
		configuration[0].S = spotPrice;
		for (int j = 1; j < n; ++j) {
			double X1 = (double)(rand() % 2 * 2 - 1);
			configuration[j].S = configuration[j - 1].S * exp(volatility * X1 * sqrtDeltaT);
		}

		double V = strike - configuration[n - 1].S;
		if (V < 0.0) V = 0.0;

		ret += V;
	}
	ret /= (double)m;
	ret *= exp(-riskFreeInterestRate * yearsToExpiry);

	return ret;
}
