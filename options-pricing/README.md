# Pricing options using Monte Carlo

Calculate price of European Vanilla and Exotic options using 
the Malliavin method explained in the 
following two papers :

#### Fournié, E. and Lasry,  J.-M. and Lebuchoux, J. and Lions, P.-L. and Touzi, N, Applications of Malliavin calculus to Monte Carlo methods in finance.

#### Fournié, E. and Lasry,  J.-M. and Lebuchoux, J. and Lions, P.-L. and Touzi, N, Applications of Malliavin calculus to Monte Carlo methods in finance. II

## Definitions

Let K be the strike price.

Let t be the time.

Let r be the risk-free interest rate.

Let sigma be the underlying volatility.

Let q be the dividend yield.

Let X be a some random variable.

Let St be the spot price.

Let Ct be the call option price.

## include folder

Contains :

    - OptionsPricingModel.h : base class
    for any model or method.
    - random_singelton.h : Practical_Class_M2MO


## Black Scholes
The folder black-scholes contains
functions for computing the Black-Scholes formulas
of European Vanilla calls and puts
(price + Greeks).

## Monte Carlo method

Based on the Malliavin calculus method explained 
in the papers mentioned above.
Two functions :

    - MalliavinEuropeanVanilla
    - MalliavinAsianExotic

which compute the price and greeks for 
any vanilla european product + an Asian exotic option.

## Main files

We saw fit to define two main files:

    - main_rd.cpp : generates simulations 
    a European vanilla product and saves 
    them in .txt files
    - main_asian.cpp : same as 
    main_rd.cpp but for Asian exotic 
    options

