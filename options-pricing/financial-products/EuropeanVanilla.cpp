//
// Created by Mohamed Taha Bennani on 06/02/2022.
//

#include "EuropeanVanilla.h"

double EuropeanVanilla::CallPayoff(double spotPrice, double strike) {
    return (spotPrice>=strike)? spotPrice-strike:0.0 ;
}

double EuropeanVanilla::DigitalPayoff(double spotPrice, double strike) {
    return (spotPrice>=strike)? 1.0:0.0 ;
}

double EuropeanVanilla::PutPayoff(double spotPrice, double strike) {
    return (strike>=spotPrice)? strike-spotPrice:0.0 ;
}