# Pricing options using Monte Carlo

Calculate price, implied volatility of European options with Black Scholes' model, Binomial model and Monte Carlo model.

## Definitions

Let K be the strike price.

Let t be the time.

Let r be the risk-free interest rate.

Let sigma be the underlying volatility.

Let X be a some random variable.

Let S(t, X) be the spot price.

Let C(S, t) be the call option price.

Let delta be round(C)/round(S), theta be round(C)/round(t) and gamma be round^2(C)/round(S)^2


Binomial Model and Monte Carlo Model
=====

Based on delta hedging and that X follows geometric Brownian motion, using different options pricing model,
we can compute the same results as Black Scholes' model.
Binomial model computes option prices with a probability tree and Monte Carlo model computes option prices by simulation. In different problem settings, like calculating American option price, where analytical solutions do not exist, these models are particularly useful.

Implied Volatility
=====
In reality, the unknown is not the option price but the underlying volatility. We can still calculate this unknown implied by the model with all sorts of numerical methods.
