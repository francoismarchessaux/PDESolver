#include <cmath>
#include "BlackScholes.h"

// Default Constructor
BlackScholes::BlackScholes() : m_option(), m_S(100.0), m_r(0.02), m_vol(0.2) { }

// Initializer Constructor
BlackScholes::BlackScholes(const Option& option, double S, double r, double vol) : m_option(option), m_S(S), m_r(r), m_vol(vol) { }

// Distribution of a standard normal distribution
double BlackScholes::normDist(double x) const { return 0.5 * (1.0 + erf(x / sqrt(2.0))); }

// Function to compute the price of a vanilla option
double BlackScholes::price() const
{
    double K = m_option.getStrike();
    double T = m_option.getExpiry();

    double numerator = log(m_S / K) + T * (m_r + 0.5 * pow(m_vol, 2));
    double denominator = m_vol * sqrt(T);

    double d1 = numerator / denominator;
    double d2 = d1 - m_vol * sqrt(T);

    if(m_option.getOptionType() == OptionType::Call)
        return m_S * normDist(d1) - K * exp(-m_r * T) * normDist(d2);
    else
        return -m_S * normDist(-d1) + K * exp(-m_r * T) * normDist(-d2);
}