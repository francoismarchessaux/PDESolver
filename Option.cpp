#include <algorithm>
#include <iostream>
#include "Option.h"

// Default Constructor
Option::Option() : m_T(1.0), m_K(110.0), m_div(0.0), m_repo(0.0), m_Type(OptionType::Call) { }

// Initializer Constructor
Option::Option(double expiry, double strike, double div, double repo, OptionType optionType) : m_T(expiry), m_K(strike), m_div(div), m_repo(repo), m_Type(optionType) { }

// Operators
std::ostream& operator<<(std::ostream& os, const Option& option)
{
    os << (option.m_Type == OptionType::Call ? "Call" : "Put") << " Option with Maturity " << option.m_T << "Y and Strike $" << option.m_K;
    return os;
}

// Accessors
OptionType Option::getOptionType() const { return m_Type; }
double Option::getStrike() const { return m_K; }
double Option::getExpiry() const { return m_T; }
double Option::getDiv() const { return m_div; }
double Option::getRepo() const { return m_repo; }

// Methods
double Option::getPayoff(double S) const
{
    if(m_Type == OptionType::Call)
        return std::max(S - m_K, 0.0);
    else
        return std::max(m_K - S, 0.0);
}

