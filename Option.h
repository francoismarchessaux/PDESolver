#pragma once
#include <iostream>

// Class to represent the type of an option
enum class OptionType {Call, Put};

// Class to represent an option with its parameters
class Option
{
    public:
        // Default Constructor
        Option();

        // Initializer Constructor
        Option(double expiry, double strike, OptionType optionType);

        // Destructor
        virtual ~Option() {};

        // Operators
        friend std::ostream& operator<<(std::ostream& os, const Option& option);

        // Accessors
        OptionType getOptionType() const;
        double getExpiry() const;
        double getStrike() const;

        // Methods
        double getPayoff(double S) const;

    private:
        OptionType m_Type;
        double m_T;
        double m_K;
};