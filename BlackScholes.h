#pragma once
#include "Option.h"

class BlackScholes
{
    public:
        // Default constructor
        BlackScholes();
        
        // Initializer Constructor
        BlackScholes(const Option& option, double S, double r, double vol);

        // Destructor
        virtual ~BlackScholes() {};

        // Methods
        double price() const;
        double normDist(double x) const;

    private:
        Option m_option;
        double m_S;
        double m_r;
        double m_vol;
};