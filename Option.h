#pragma once
#include <iostream>

// Class to represent the type of an option
enum OptionType {Call, Put};

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
        OptionType getOptionType();
        double getExpiry();
        double getStrike();

        // Methods
        double getPayoff(double S);

    private:
        OptionType m_Type;
        double m_T;
        double m_K;
};