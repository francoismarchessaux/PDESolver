#include <algorithm>
#include <iostream>
#include "Option.h"

// Default Constructor
Option::Option() : m_T(1.0), m_K(110.0), m_Type(OptionType::Call) {}

// Initializer Constructor
Option::Option(double expiry, double strike, OptionType optionType) : m_T(expiry), m_K(strike), m_Type(optionType) { }

// Operators
std::ostream& operator<<(std::ostream& os, const Option& option)
{
    os << (option.m_Type == OptionType::Call ? "Call" : "Put") << " Option with Maturity " << option.m_T << "Y and Strike $" << option.m_K;
    return os;
}

// Accessors
OptionType Option::getOptionType() { return this->m_Type; }
double Option::getExpiry() { return this->m_T; }
double Option::getStrike() { return this->m_K; }

// Methods
double Option::getPayoff(double S)
{
    if(this->m_Type == OptionType::Call)
        return std::max(S - this->m_K, 0.0);
    else
        return std::max(this->m_K - S, 0.0);
}

