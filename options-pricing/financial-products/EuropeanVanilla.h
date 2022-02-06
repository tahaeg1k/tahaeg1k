//
// Created by Mohamed Taha Bennani on 06/02/2022.
//

#ifndef PRICING_EUROPEANVANILLA_H
#define PRICING_EUROPEANVANILLA_H

#endif //PRICING_EUROPEANVANILLA_H

class EuropeanVanilla {
public :
    double CallPayoff(double spotPrice, double strike);
    double DigitalPayoff(double spotPrice, double strike);
    double PutPayoff(double spotPrice, double strike);
};
