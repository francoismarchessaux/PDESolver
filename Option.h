#pragma once
#include <iostream>

enum class OptionType {Call, Put};

class Option
{
    public:
        // Default Constructor
        Option();

        // Initializer Constructor
        Option(double expiry, double strike, double div, double repo, OptionType optionType);

        // Destructor
        ~Option() {};

        // Operators
        friend std::ostream& operator<<(std::ostream& os, const Option& option);

        // Accessors
        OptionType getOptionType() const;
        double getExpiry() const;
        double getStrike() const;
        double getDiv() const;
        double getRepo() const;

        // Methods
        double getPayoff(double S) const;

    private:
        OptionType m_Type;
        double m_T;
        double m_K;
        double m_div;
        double m_repo;
};