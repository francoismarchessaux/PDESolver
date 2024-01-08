#include <cmath>
#include "BlackScholes.h"
#include "Option.h"

// Default Constructor
BlackScholes::BlackScholes() : m_option(), m_S(100.0), m_r(0.02), m_vol(0.2) { }

// Initializer Constructor
BlackScholes::BlackScholes(const Option& option, double S, double r, double vol) : m_option(option), m_S(S), m_r(r), m_vol(vol) { }

// Distribution of a standard normal distribution
double BlackScholes::normDist(double x) const { return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0))); }

// Function to compute the price of a vanilla option
double BlackScholes::price() const
{
    double K = m_option.getStrike();
    double T = m_option.getExpiry();
    double div = m_option.getDiv();
    double repo = m_option.getRepo();

    double numerator = std::log(m_S / K) + T * (m_r - div - repo + 0.5 * std::pow(m_vol, 2));
    double denominator = m_vol * std::sqrt(T);

    double d1 = numerator / denominator;
    double d2 = d1 - m_vol * std::sqrt(T);

    if(m_option.getOptionType() == OptionType::Call)
        return m_S * normDist(d1) * std::exp(T * (-repo - div)) - K * std::exp(-m_r * T) * normDist(d2);
    else
        return -m_S * normDist(-d1)  * std::exp(T * (-repo - div)) + K * std::exp(-m_r * T) * normDist(-d2);
}