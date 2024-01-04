#include <cmath>
#include "BlackScholes.h"

// Default Constructor
BlackScholes::BlackScholes() : m_option(), m_S(100.0), m_r(0.02), m_vol(0.2) { }

// Constructor
BlackScholes::BlackScholes(Option* option, double S, double r, double vol)
{
    m_option = option;
    m_S = S;
    m_r = r;
    m_vol = vol;
}

// Distribution of a standard normal distribution
double BlackScholes::normDist(double x)
{
    return 0.5 * (1 + erf(x / sqrt(2)));
}

// Function to compute the price of a vanilla option
double BlackScholes::operator()()
{
    double K = this->m_option->getStrike();
    double T = this->m_option->getExpiry();

    double numerator = log(this->m_S / K) + T * (this->m_r + 0.5 * pow(this->m_vol, 2));
    double denominator = this->m_vol * sqrt(T);

    double d1 = numerator / denominator;
    double d2 = d1 - this->m_vol * sqrt(T);

    if(this->m_option->getOptionType() == OptionType::Call)
        return this->m_S * this->normDist(d1) - K * exp(-this->m_r * T) * this->normDist(d2);
    else
        return -m_S * this->normDist(-d1) + K * exp(-m_r * T) * this->normDist(-d2);
}