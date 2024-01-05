#pragma once
#include <functional>
#include <vector>
#include "Matrix.h"
#include "Option.h"

using namespace std;

enum PDEType {Backward, Forward};

class PDE
{
    public:
        // Default Constructor
        PDE();

        // Initializer Constructors
        PDE(size_t nTimeSteps, double minTime, double maxTime, size_t nSpaceSteps, double multiplier, double vol, double spot);
        PDE(vector<double>& timeGrid, vector<double> spaceGrid);

        // Copy Constructor
        PDE(const PDE& rhs); 

        // Destructor
        virtual ~PDE() {};

        // Assignement operator
        PDE& operator=(const PDE& rhs);

        // Accessors
        vector<double> getTimeGrid() const;
        vector<double> getSpaceGrid() const;
        Matrix getSolution() const;

        // Mutators
        void setTimeGrid(size_t nTimeSteps, double minTime, double maxTime);
        void setTimeGrid(const vector<double>& timeGrid);
        void setSpaceGrid(size_t nSpaceSteps, double multiplier, double vol, double spot);
        void setSpaceGrid(const vector<double>& spaceGrid);

        // Partial derivatives
        double partial_t(size_t i);
        double partial_x(size_t i);

        // Matrices
        Matrix computeP(size_t i, size_t m, double theta);
        Matrix computeQ(size_t i, size_t m, double theta);
        vector<double> computeV(size_t i, size_t m, double theta);

        // PDE methods
        void setPDE(function<double (double, double)>& a, function<double (double, double)>& b, function<double (double, double)>& c, function<double (double, double)>& d);
        vector<double> createTerminalCondition(Option* option);
        void setProblem(vector<double> leftBoundary, vector<double> rightBoundary, vector<double> terminalCondition);
        void resolve();
        double solution(double spot);

    private:
        vector<double> m_timeGrid;
        vector<double> m_spaceGrid;
        vector<double> m_leftBoundary;
        vector<double> m_rightBoundary;
        vector<double> m_terminalCondition;
        function<double (double, double)> m_a;
        function<double (double, double)> m_b;
        function<double (double, double)> m_c;
        function<double (double, double)> m_d;
        Matrix m_Solution;
};