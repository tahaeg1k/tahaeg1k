//
//  random-singleton.cpp
//  Practical_Class_M2MO
//
//  Created by Noufel frikha on 17/11/2020.
//  Copyright © 2020 Noufel frikha. All rights reserved.
//
#include "random_singleton.h"

// random-singleton.cpp (Singleton)
// Implementation de la classe Random / Implementation of the class Random

#include <cfloat> //DBL_EPSILON
#include <iostream>
using namespace std;

inline static double sqr(double x) {return x*x;}

//Ces constantes sont requises par l'algorithme du generateur aleatoire
//These consts are required by the algorithm
const int Random::NTAB = 32;

const long Random::IM1 = 2147483563;
const long Random::IM2 = 2147483399;
const long Random::IMM1 = Random::IM1 - 1;
const double Random::AM = 1.0/Random::IM1;
const int Random::IA1 = 40014;
const int Random::IA2 = 40692;
const int Random::IQ1 = 53668;
const int Random::IQ2 = 52774;
const int Random::IR1 = 12211;
const int Random::IR2 = 3791;
const int Random::NDIV = 1 + Random::IMM1/Random::NTAB;
const double Random::EPS = DBL_EPSILON;
const double Random::RNMX = 1.0-Random::EPS;

long  Random::idum   = 0;
long  Random::idum2  = 123456789L;
long  Random::iy     = 0;
long* Random::iv     = 0;

//definition du singleton
Random Random::Singleton(0);

void Random::Randomize(long seed)
{
    idum = (seed <= 0) ? (seed == 0 ? 1 : -seed) : seed; //be sure to prevent idum=0

    idum2 = idum;
    if (idum2 < 0) idum2 = -idum2;

    long k = 0;
    for(int j = NTAB+7 ; j >= 0 ; j--)
    {
        k = idum/IQ1;
        idum = IA1*(idum - k*IQ1) - k*IR1;
        if (idum < 0)
            idum += IM1;

        if (j < NTAB)
            iv[j] = idum;
    }

    iy = iv[0];
}
//fin Randomize()

//Cette fonction renvoie un double pseudo-aleatoire uniformement dans ]0;1[
//Il s'agit de l'algorithme de L'Ecuyer avec melange de Bays-Durham
//This function returns a double, taken uniformly in ]0;1[
//It is the algorithm of L'Ecuyer with a Bays-Durham shuffle
double Random::theRandom(void)
{
    long k = idum/IQ1;
    idum = IA1*(idum - k*IQ1) - k*IR1;
    if (idum < 0)
        idum +=IM1;

    k = idum2/IQ2;
    idum2 = IA2*(idum2 -k*IQ2) - k*IR2;
    if (idum2 < 0)
        idum2 +=IM2;

    long j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;
    if (iy < 1)
        iy += IMM1;

    double temp = AM*iy;
    if (temp >= RNMX)
        return RNMX;  //empeche de renvoyer 1 / prevents from returning 1
    else
        return temp;
}
//fin theRandom()

//Renvoie un double de la Gaussienne de moyenne et d'ecartype specifies
//Returns a double taken on a Gaussian with specified mean and standard deviation
double Random::Gaussian(double mean, double standardDeviation)
{
    const int NbTirages = 12; //augmenter pour une meilleure precision. 12 est bien.
    //increase for better precision. 12 works fine.
    double valeur = 0;
    for(int i=0 ; i < NbTirages ; ++i)
        valeur += Random::Uniform();

    //on recentre la somme / centering the sum
    valeur -= NbTirages/2;

    //on etale suivant l'ecartype / spread with standard deviation
    //le 12 n'a rien a voir avec NbTirages, mais explique pourquoi justement, on prend souvent
    //NbTirages = 12
    //the 12 is not related to NbTirages, but it explains why it is often chosen that NbTirages=12
    valeur *= (NbTirages == 12) ? standardDeviation
                                : sqrt(12/static_cast<double>(NbTirages))*standardDeviation;

    //on centre sur la moyenne / debias
    valeur += mean;

    return valeur;
}

//fin Gaussian()

double Random::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/sqrt(2))/2;
}

double Random::normalPDF(double x) // n(x)
{
    return std::exp(-x*x/2)/ std::sqrt(2*M_PI);
}

double Random::inverse_normalCDF(double x)
{
    double c0 = 2.515517;
    double c1 = 0.802853;
    double c2 = 0.010328;
    double d1 = 1.432788;
    double d2 = 0.189269;
    double d3 = 0.001308;
    double signe, t;

    if (x>0.5)
    {
        signe = +1.0; x=1.0-x;
    }
    else signe = -1.0;

    t = sqrt(-2.0 * log(x));
    double result = signe * (t-((c2*t+c1)*t+c0)/(1.0+t*(d1+t*(d2+d3*t))));

    return result;
}

double Random::Payoff_Call(double S, double K)
{
    return (S >= K) ? S - K : 0.0;
}

double Random::Payoff_Digital_Call(double S, double K)
{
    return (S >= K) ? 1 : 0.0;
}

double Random::Payoff_Best_of_Call(double ST1, double ST2, double K)
{
    return (max(ST1,ST2)>=K) ? max(ST1,ST2) - K : 0.0;
}

double Random::Call_price_BS(double S, double r, double T, double K, double sigma)
{
    double d_1 = 1 / (sigma * sqrt(T)) * ((log(S) - log(K) + (r + pow(sigma, 2)/2) * T));
    return S * normalCDF(d_1) - K * exp(- r * T) * normalCDF(d_1 - sigma * sqrt(T));
}

int Random::Binomial(double prob)
{
    int result = 0;
    double unif = Random::Uniform();
    if( unif < prob)
    {
        result = 1;
    }
    return result;
}

int Random::Bernoulli(int N, double prob)
{
    int result = 0;
    for(int n=0; n<N; n++)
    {
        result += Random::Binomial(prob);
    }
    return result;
}

double Random::Exponentielle(double lambda)
{
    double result;
    double unif = Random::Uniform();
    result = (-1/lambda) * log(unif);
    return result;
}

double Random::Cauchy(double x0, double lambda)
{
    double result;
    double unif = Random::Uniform();
    result = lambda * tan(M_PI*(unif-0.5)) + x0;
    return result;
}

double * Random::BoxMuller(double m, double sigma)
{
    double * result = new double [2];
    double u = Random::Uniform();
    double expo = Random::Exponentielle(0.5);

    result[0] = m + sigma * sqrt(expo) * cos(2 * M_PI * u);

    result[1] = m + sigma * sqrt(expo) * sin(2 * M_PI * u);

    return result;
}
double Random::Gamma(double a)
{
    // Méthode d'acceptation - rejet
    double result;
    double unif1;
    double unif2;
    double Y;
    double q_Y;
    double seuil  = exp(1)/(a + exp(1));

    int arret = 0;

    while(arret == 0)
    {
        unif1 = Random::Uniform();
        unif2 = Random::Uniform();

        if(unif2 < seuil)
            Y = pow(unif2/seuil,1/a);
        else
            Y = -log((1-unif2)*(1/(seuil*exp(1))));
        if(Y< 1)
            q_Y = exp(-Y);
        else
            q_Y = pow(Y, a-1);

        if(unif1 < q_Y)
            arret = 1;
    }
    result = Y;
    return result;
}
double * Random::Inverse_Gaussian(double lambda, double * mu, int length_mu)
{
    double z= Random::Gaussian(0,1);
    double u=  Random::Uniform(0,1);
    double y = z*z;
    double * x = new double [length_mu];
    for(int i=0; i<length_mu; ++i){
        x[i]=mu[i] + (0.5*mu[i]/lambda)*(y*mu[i] +sqrt(4*lambda*mu[i]*y+mu[i]*mu[i]*y*y));
        if (u > mu[i]/(mu[i]+x[i])) x[i]=mu[i]*mu[i]/x[i];
    }
    return x;
}

double * Random::Euler_scheme_BS(double s0, double sigma, double mu, double T, int M, int N)
{
    double * result = new double[2];

    double scur;
    double ST;
    double delta = T/(double(N));
    double G;
    double increm_brown;
    double brownien_cur;
    double error_local_l1 = 0;
    double error_local_l2 = 0;

    for(int k=0; k< M; k++)
    {
        scur = s0;
        brownien_cur = 0;
        for(int i=0; i<N; i++)
        {
            G = Random::Gaussian();
            increm_brown = sqrt(delta) * G;
            brownien_cur += increm_brown;
            scur += mu * scur * delta + sigma * scur * increm_brown;
        }
        ST = s0 * exp((mu-0.5*sigma*sigma)*T + sigma* brownien_cur);
        error_local_l1 += abs(ST-scur);
        error_local_l2 += (ST-scur) * (ST-scur);
    }
    error_local_l1 /= double(M);
    error_local_l2 /= double(M);

    error_local_l1 *= sqrt(double(N));
    error_local_l2 *= double(N);

    result[0] = error_local_l1;
    result[1] = error_local_l2;
    return result;
}
double Random::Weak_error_BS(double s0, double sigma, double r, double T, double K, int M, int N)
{
    double result;
    double price_MC = 0;
    double price_BS = Random::Call_price_BS(s0, r, T, K, sigma);
    double scur;
    double delta = T/(double(N));
    double G;
    double increm_brown;

    for(int k=0; k< M; k++)
    {
        scur = s0;
        for(int i=0; i<N; i++)
        {
            G = Random::Gaussian();
            increm_brown = sqrt(delta) * G;
            scur += r * scur * delta + sigma * scur * increm_brown;
        }
        price_MC += Random::Payoff_Call(scur, K);
    }
    price_MC /= double(M);
    price_MC *= exp(-r*T);
    result = abs(price_BS-price_MC) * double(N);

    return result;
}


double * Random::Milstein_scheme_BS(double s0, double sigma, double mu, double T, int M, int N)
{
    double * result = new double[2];

    double scur;
    double ST;
    double delta = T/(double(N));
    double G;
    double increm_brown;
    double brownien_cur;
    double error_local_l1 = 0;
    double error_local_l2 = 0;

    for(int k=0; k< M; k++)
    {
        scur = s0;
        brownien_cur = 0;
        for(int i=0; i<N; i++)
        {
            G = Random::Gaussian();
            increm_brown = sqrt(delta) * G;
            brownien_cur += increm_brown;
            scur += mu * scur * delta + sigma * scur * increm_brown + 0.5 * sigma * sigma * scur * (increm_brown*increm_brown-delta);
        }
        ST = s0 * exp((mu-0.5*sigma*sigma)*T + sigma* brownien_cur);
        error_local_l1 += abs(ST-scur);
        error_local_l2 += (ST-scur) * (ST-scur);
    }
    error_local_l1 /= double(M);
    error_local_l2 /= double(M);

    error_local_l1 *= double(N);
    error_local_l2 *= double(N) * double(N);

    result[0] = error_local_l1;
    result[1] = error_local_l2;
    return result;
}

double Random::VaR_BS(double alpha, int M, double theta_init, double S, double r, double T, double sigma)
{
    double VaR_theorique = S*exp((r-0.5*sigma*sigma)*T + sigma * sqrt(T) * inverse_normalCDF(alpha));

    double theta_cur = theta_init;
    double G, ST, H;

    for(int i=0; i< M; i++)
    {
        G = Random::Gaussian();
        ST = S*exp((r-0.5*sigma*sigma)*T + sigma * sqrt(T) * G);
        H = (ST <= theta_cur) ? 1-alpha : -alpha;

        theta_cur -= (1/(10+pow(i+1, 3/4))) * H;
        cout << theta_cur << " ; " << "VaR theorique = " << VaR_theorique << endl;
    }
    return theta_cur;
}
double * Random::IS_alg_sto(double theta0, int M, double S0, double L, double sigma, double r, double T)
{
    double theta_cur = theta0;
    double G, ST, payoff, payoff_IS;
    double * result = new double[3];
    double prix_MC = 0;
    double new_E2 = 0;
    double payoff_theta;
    double prix_MC_sans_IS = 0;
    double old_E2 = 0;

    for(int i=0; i<M; i++)
    {
        G = Random::Gaussian();
        ST = S0*exp((r-0.5*sigma*sigma)*T + sigma *sqrt(T)*(G-theta_cur));
        payoff_IS = (ST >= L) ? 1 : 0;


        ST = S0*exp((r-0.5*sigma*sigma)*T + sigma *sqrt(T)*(G+theta_cur));
        payoff_theta = (ST >= L) ? 1 : 0;

        prix_MC += payoff_theta * exp(-theta_cur*G - 0.5 * theta_cur*theta_cur);
        new_E2 += payoff_theta * payoff_theta* exp(-2*theta_cur*G - theta_cur*theta_cur);

        ST = S0*exp((r-0.5*sigma*sigma)*T + sigma *sqrt(T)*G);
        payoff = (ST >= L) ? 1 : 0;

        prix_MC_sans_IS += payoff;
        old_E2 += payoff * payoff;

        theta_cur -= (1/(10+pow(i+1, 3/4))) * payoff_IS*payoff_IS*(2*theta_cur-G);
        cout << " theta = " << theta_cur << endl;


    }

    prix_MC /= double(M);
    new_E2 /= double(M);
    prix_MC_sans_IS /= double(M);
    old_E2 /= double(M);

    double variance_IS = new_E2 - prix_MC*prix_MC;
    double variance_sans_IS = old_E2 - prix_MC_sans_IS *prix_MC_sans_IS;

    cout << " prix MC sans IS = " << prix_MC_sans_IS << " prix MC avec IS = " << prix_MC << endl;
    cout << " variance sans IS = " << variance_sans_IS << " variance avec IS = " << variance_IS << endl;
    result[0] = theta_cur;
    result[1] = prix_MC;
    result[2] = prix_MC_sans_IS;
    return result;
}

double * Random::Newton_Raphson(double *G, double theta0, double epsilon, int M, int length_G, double S0, double K, double sigma, double r, double T)
{
    double theta_cur = theta0;
    double grad_comp = 0;
    double hess_comp = 0;
    double payoff_carre;
    double ST;
    double ST_translate;
    double ST_inv;
    bool prec_achieved = true;
    double EM_IS = 0;
    double EM_standard =0;
    double EM2_IS =0;
    double EM2_standard =0;
    double EM_AV = 0;
    double EM2_AV = 0;

    double * result = new double[5];

    while(prec_achieved)
    {

        for(int i=0; i< M; i++)
        {
            ST = S0 * exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)*G[i]);
            payoff_carre = Payoff_Digital_Call(ST, K) * Payoff_Digital_Call(ST, K);
            grad_comp += (theta_cur-G[i]) * payoff_carre * exp(-theta_cur*G[i]+ 0.5*theta_cur*theta_cur);
            hess_comp += (1+(theta_cur-G[i])*(theta_cur-G[i]))* payoff_carre * exp(-theta_cur*G[i]+ 0.5*theta_cur*theta_cur);
        }
        grad_comp /= double(M);
        hess_comp /= double(M);
        theta_cur -= (1/hess_comp)*grad_comp;
        if(abs(grad_comp) < epsilon)
        {
            prec_achieved = false;
        }
        grad_comp = 0;
        hess_comp = 0;

    }
    for(int i=0; i< M; i++)
    {
        ST = S0 * exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)* G[i]);
        ST_translate = S0 * exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)*(G[i]+theta_cur));
        EM_IS += Payoff_Digital_Call(ST_translate, K) *  exp(-theta_cur*G[i] - 0.5*theta_cur*theta_cur);
        EM2_IS += Payoff_Digital_Call(ST_translate, K) * Payoff_Digital_Call(ST_translate, K)* exp(-2*theta_cur*G[i] - theta_cur*theta_cur);
        EM_standard += Payoff_Digital_Call(ST, K);
        EM2_standard += Payoff_Digital_Call(ST, K) * Payoff_Digital_Call(ST, K);
        ST_inv = S0 * exp((r-0.5*sigma*sigma)*T - sigma*sqrt(T)* G[i]);
        EM_AV +=  0.5* (Payoff_Digital_Call(ST, K) +  Payoff_Digital_Call(ST_inv, K));
        EM2_AV += 0.25*(Payoff_Digital_Call(ST, K) +  Payoff_Digital_Call(ST_inv, K))*(Payoff_Digital_Call(ST, K) +  Payoff_Digital_Call(ST_inv, K));
    }
    EM_IS /= double(M);
    EM_standard /= double(M);
    EM2_IS /= double(M);
    EM2_standard /= double(M);
    EM_AV /= double(M);
    EM2_AV /= double(M);

    double variance_IS = EM2_IS-EM_IS*EM_IS;
    double variance_standard  = EM2_standard -EM_standard*EM_standard;
    double variance_AV = EM2_AV - EM_AV*EM_AV;

    result[0] = theta_cur;
    result[1] = EM_IS;
    result[2] = 1.96*sqrt(variance_IS)/sqrt(double(M));
    result[3] = EM_AV;
    result[4] = variance_standard/variance_AV;
    result[5] = variance_standard/variance_IS;
    return result;
}
//fin Newton_Raphson_1
double * Random::Newton_Raphson_Best_of_Call(double *G1, double *G2, double theta01, double theta02, double epsilon, int M, int length_G, double S01, double S02, double K, double sigma1, double sigma2, double rho, double r, double T)
{
    double theta_cur1 = theta01;
    double theta_cur2 = theta02;
    double grad_comp1 = 0;
    double grad_comp2 = 0;
    double hess_comp11, hess_comp12, hess_comp21, hess_comp22 = 0;
    double det = 0;
    double inv_hess_time_grad1, inv_hess_time_grad2 = 0;
    double payoff_carre;
    double ST1;
    double ST1_translate;
    double ST2;
    double ST2_translate;
    double ST1_inv, ST2_inv;
    bool prec_achieved = true;
    double EM_IS = 0;
    double EM_standard =0;
    double EM2_IS =0;
    double EM2_standard =0;
    double EM_AV = 0;
    double EM2_AV = 0;

    double * result = new double[5];

    while(prec_achieved)
    {

        for(int i=0; i< M; i++)
        {
            ST1 = S01 * exp((r-0.5*sigma1*sigma1)*T + sigma1*sqrt(T)*G1[i]);

            ST2 = S02 * exp((r-0.5*sigma2*sigma2)*T + sigma2*sqrt(T)*(rho*G1[i]+sqrt(1-rho*rho)*G2[i]));

            payoff_carre = Payoff_Best_of_Call(ST1, ST2, K) * Payoff_Best_of_Call(ST1, ST2, K);

            grad_comp1 += (theta_cur1-G1[i]) * payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

            grad_comp2 += (theta_cur2-G2[i]) * payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

            hess_comp11 += (1+(theta_cur1-G1[i])*(theta_cur1-G1[i]))* payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

            hess_comp12 += ((theta_cur1-G1[i])*(theta_cur2-G2[i]))* payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

            hess_comp21 += ((theta_cur2-G2[i])*(theta_cur1-G1[i]))* payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));


            hess_comp22 += (1+(theta_cur2-G2[i])*(theta_cur2-G2[i]))* payoff_carre * exp(-theta_cur1*G1[i]-theta_cur2*G2[i]+ 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

        }
        grad_comp1 /= double(M);
        grad_comp2 /= double(M);
        hess_comp11 /= double(M);
        hess_comp12 /= double(M);
        hess_comp21 /= double(M);
        hess_comp22 /= double(M);

        det = hess_comp11 * hess_comp22 - hess_comp12 * hess_comp21;

        inv_hess_time_grad1 = (1/det)* (hess_comp22*grad_comp1-hess_comp12*grad_comp2);

        inv_hess_time_grad2 = (1/det)* (hess_comp11*grad_comp2-hess_comp21*grad_comp1);


        theta_cur1 -= inv_hess_time_grad1;

        theta_cur2 -= inv_hess_time_grad2;

        cout << "theta1 = " << theta_cur1 << endl;

        cout << "theta2 = " << theta_cur2 << endl;

        if(sqrt((grad_comp1*grad_comp1) + (grad_comp2*grad_comp2)) < epsilon)
        {
            prec_achieved = false;
        }

        grad_comp1  = 0;
        grad_comp2  = 0;
        hess_comp11 = 0;
        hess_comp12 = 0;
        hess_comp21 = 0;
        hess_comp22 = 0;

    }
    for(int i=0; i< M; i++)
    {
        ST1 = S01 * exp((r-0.5*sigma1*sigma1)*T + sigma1*sqrt(T)*G1[i]);

        ST2 = S02 * exp((r-0.5*sigma2*sigma2)*T + sigma2*sqrt(T)*(rho*G1[i]+sqrt(1-rho*rho)*G2[i]));

        ST1_translate = S01 * exp((r-0.5*sigma1*sigma1)*T + sigma1*sqrt(T)*(G1[i]+theta_cur1));

        ST2_translate = S02 * exp((r-0.5*sigma2*sigma2)*T + sigma2*sqrt(T)*(rho*(G1[i]+theta_cur1)+sqrt(1-rho*rho)*(G2[i]+theta_cur2)));
        ;

        EM_IS += Payoff_Best_of_Call(ST1_translate, ST2_translate, K)  *  exp(-theta_cur1*G1[i]-theta_cur2*G2[i] - 0.5*(theta_cur1*theta_cur1+theta_cur2*theta_cur2));

        EM2_IS += Payoff_Best_of_Call(ST1_translate, ST2_translate, K) * Payoff_Best_of_Call(ST1_translate, ST2_translate, K) *  exp(-2*theta_cur1*G1[i]-2*theta_cur2*G2[i] - (theta_cur1*theta_cur1+theta_cur2*theta_cur2));

        EM_standard += Payoff_Best_of_Call(ST1, ST2, K)  ;
        EM2_standard += Payoff_Best_of_Call(ST1, ST2, K)* Payoff_Best_of_Call(ST1, ST2, K);
        ST1_inv = S01 * exp((r-0.5*sigma1*sigma1)*T - sigma1*sqrt(T)* G1[i]);
        ST2_inv = S02 * exp((r-0.5*sigma2*sigma2)*T - sigma2*sqrt(T)* (rho*G1[i]+sqrt(1-rho*rho)*G2[i]));
        EM_AV +=  0.5* (Payoff_Best_of_Call(ST1, ST2, K) +  Payoff_Best_of_Call(ST1_inv, ST2_inv, K));
        EM2_AV += 0.25*(Payoff_Best_of_Call(ST1, ST2, K) +  Payoff_Best_of_Call(ST1_inv, ST2_inv, K))*(Payoff_Best_of_Call(ST1, ST2, K) +  Payoff_Best_of_Call(ST1_inv, ST2_inv, K));
    }
    EM_IS /= double(M);
    EM_standard /= double(M);
    EM2_IS /= double(M);
    EM2_standard /= double(M);
    EM_AV /= double(M);
    EM2_AV /= double(M);

    double variance_IS = EM2_IS-EM_IS*EM_IS;
    double variance_standard  = EM2_standard -EM_standard*EM_standard;
    double variance_AV = EM2_AV - EM_AV*EM_AV;

    result[0] = theta_cur1;
    result[1] = EM_IS;
    result[2] = 1.96*sqrt(variance_IS)/sqrt(double(M));
    result[3] = EM_AV;
    result[4] = variance_standard/variance_AV;
    result[5] = variance_standard/variance_IS;
    return result;
}

double * Random::Normal_Inverse_Gaussian(double alpha, double beta, double delta, double mu, double * shift, int length_shift)
{
    double z = Random::Gaussian(0,1);
    double * gamma= new double[length_shift];
    double * mu_IG= new double[length_shift];
    double * res= new double[length_shift];
    //Convention: le shift se fait en soustrayant!
    for(int i=0; i< length_shift; ++i){
        gamma[i]=sqrt(alpha*alpha-(beta-shift[i])*(beta-shift[i]));
        mu_IG[i]=delta/gamma[i];
    }
    double * y =Random::Inverse_Gaussian(delta*delta,mu_IG, length_shift);
    for(int i=0; i< length_shift; ++i) res[i]=mu + (beta-shift[i])*y[i] + sqrt(y[i])*z;

    //delete gamma, mu_IG; // warning from GCC
    delete[] gamma;
    delete[] mu_IG;
    return res;
}

//fin Normal_Inverse_Gaussian()

//fin random-singleton.cpp