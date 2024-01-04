#pragma once
#include "Option.h"

class BlackScholes
{
    public:
        // Default constructor
        BlackScholes();
        
        // Constructor
        BlackScholes(Option* option, double S, double r, double vol);

        // Destructor
        virtual ~BlackScholes() {};

        // Methods
        double operator()();
        double normDist(double x);

    private:
        Option* m_option;
        double m_S;
        double m_r;
        double m_vol;
};