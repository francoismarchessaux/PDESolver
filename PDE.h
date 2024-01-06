#pragma once
#include <functional>
#include "Matrix.h"
#include "Option.h"

class PDE
{
    public:
        // Initializer Constructors
        PDE(Option option, size_t nTimeSteps, double T, size_t nSpaceSteps, double multiplier, double vol, double spot, double r, std::function<double (double, double)>& a, std::function<double (double, double)>& b, std::function<double (double, double)>& c, std::function<double (double, double)>& d);

        // Destructor
        virtual ~PDE() {};

        // Accessors
        std::vector<double> getTimeGrid() const;
        std::vector<double> getSpaceGrid() const;
        
        // Mutators
        void setTimeGrid(size_t nTimeSteps, double T);
        void setTimeGrid(const std::vector<double>& timeGrid);
        void setSpaceGrid(size_t nSpaceSteps, double multiplier, double vol, double spot);
        void setSpaceGrid(const std::vector<double>& spaceGrid);
        void setTerminalCondition();
        void setBoundaries(size_t nTimeSteps, double r);

        // Partial derivatives
        double partial_t(size_t i) const;
        double partial_x(size_t i) const;

        // Matrices computation
        Matrix computeP(size_t i, size_t m, double theta);
        Matrix computeQ(size_t i, size_t m, double theta);
        std::vector<double> computeV(size_t i, size_t m, double theta);

        // PDE methods
        void resolve();
        double solution(double spot);

    private:
        Option m_Option;
        std::vector<double> m_timeGrid;
        std::vector<double> m_spaceGrid;
        std::vector<double> m_leftBoundary;
        std::vector<double> m_rightBoundary;
        std::vector<double> m_terminalCondition;
        std::function<double (double, double)> m_a;
        std::function<double (double, double)> m_b;
        std::function<double (double, double)> m_c;
        std::function<double (double, double)> m_d;
        Matrix m_Solution;
};