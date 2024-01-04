#include <vector>
#include <functional>
#include "PDE.h"
#include "Matrix.h"

using namespace std;

// Overriding + operator for vectors
vector<double> operator+(const vector<double>& v1, const vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vector addition: vectors must have the same size");
    }

    vector<double> res(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        res[i] = v1[i] + v2[i];
    }
    return res;
}

// Default Constructor
PDE::PDE()
{
    // Initialize the time grid
    this->setTimeGrid(50, 0.0, 1.0);

    // Initialize the space grid
    this->setSpaceGrid(100, 2.0, 0.5, 100.0);

    // Initialize solution matrix
    this->m_Solution = Matrix(50, 100);
}

// Constructors
PDE::PDE(size_t nTimeSteps, double minTime, double maxTime, size_t nSpaceSteps, double multiplier, double vol, double strike)
{
    // Initialize the time grid
    this->setTimeGrid(nTimeSteps, minTime, maxTime);

    // Initialize the space grid
    this->setSpaceGrid(nSpaceSteps, multiplier, vol, strike);

    // Initialize solution matrix
    this->m_Solution = Matrix(nTimeSteps - 1, nSpaceSteps - 1);
}

PDE::PDE(vector<double>& timeGrid, vector<double> spaceGrid)
{
    // Initialize the time grid
    this->m_timeGrid = timeGrid;
    sort(this->m_timeGrid.begin(), this->m_timeGrid.end());

    // Initialize the space grid
    this->m_spaceGrid = spaceGrid;
    sort(this->m_spaceGrid.begin(), this->m_spaceGrid.end());

    // Initialize solution matrix
    this->m_Solution = Matrix(timeGrid.size() - 1, spaceGrid.size() - 1);
}

// Copy Constructor
PDE::PDE(const PDE& rhs)
{
    this->m_timeGrid = rhs.m_timeGrid;
    this->m_spaceGrid = rhs.m_spaceGrid;
    this->m_leftBoundary = rhs.m_leftBoundary;
    this->m_rightBoundary = rhs.m_rightBoundary;
    this->m_terminalCondition = rhs.m_terminalCondition;
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
    this->m_c = rhs.m_c;
    this->m_d = rhs.m_d; 
    this->m_Solution = rhs.m_Solution;
}

// Assignement operator
PDE& PDE::operator=(const PDE& rhs)
{
    if (this == &rhs)
        return *this;

    this->m_timeGrid = rhs.m_timeGrid;
    this->m_spaceGrid = rhs.m_spaceGrid;
    this->m_leftBoundary = rhs.m_leftBoundary;
    this->m_rightBoundary = rhs.m_rightBoundary;
    this->m_terminalCondition = rhs.m_terminalCondition;
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
    this->m_c = rhs.m_c;
    this->m_d = rhs.m_d; 
    this->m_Solution = rhs.m_Solution;

    return *this;
}

// Accessors
vector<double> PDE::getTimeGrid() const { return this->m_timeGrid; }
vector<double> PDE::getSpaceGrid() const { return this->m_spaceGrid; }
Matrix PDE::getSolution() const { return this->m_Solution; }

// Mutators
void PDE::setTimeGrid(size_t nTimeSteps, double minTime, double maxTime)
{
    double deltaTime = (maxTime - minTime) / static_cast<double>(nTimeSteps);

    for(size_t i = 0; i <= nTimeSteps; i++)
        this->m_timeGrid.push_back(minTime + i * deltaTime);
}

void PDE::setTimeGrid(const vector<double>& timeGrid)
{
    this->m_timeGrid = timeGrid;
    sort(this->m_timeGrid.begin(), this->m_timeGrid.end());
}

void PDE::setSpaceGrid(size_t nSpaceSteps, double multiplier, double vol, double strike)
{
    double lowerBoundary = strike - multiplier * vol * strike;
    double upperBoundary = strike + multiplier * vol * strike;
    double deltaSpace = (upperBoundary - lowerBoundary) / static_cast<double>(nSpaceSteps);

    for(size_t i = 0; i <= nSpaceSteps; i++)
        this->m_spaceGrid.push_back(lowerBoundary + i * deltaSpace);
}

void PDE::setSpaceGrid(const vector<double>& spaceGrid)
{
    this->m_spaceGrid = spaceGrid;
    sort(this->m_spaceGrid.begin(), this->m_spaceGrid.end());
}

// Partial derivatives
double PDE::partial_t(size_t i) { return this->m_timeGrid[i + 1] - this->m_timeGrid[i]; }
double PDE::partial_x(size_t i) { return this->m_spaceGrid[i + 1] - this->m_spaceGrid[i]; }

// Matrices
Matrix PDE::computeP(size_t i, size_t m, double theta)
{
    Matrix P(m - 1, m - 1);
    double t_i = this->m_timeGrid[i];

    for(size_t j = 0; j < m - 1; j++)
    {
        double x_i = this->m_spaceGrid[j];
        double x_ip1 = this->m_spaceGrid[j+1];
        P[j][j] = this->m_a(t_i, x_i) - ((1.0 / this->partial_t(i)) + ((2.0 * theta * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2)));
        
        if(j + 1 < m - 1)
        {
            P[j][j + 1] = (this->m_b(t_i, x_i)/(2 * this->partial_x(i))) + ((theta * this->m_c(t_i, x_i))/pow(this->partial_x(i), 2));
            P[j + 1][j] = -(this->m_b(t_i, x_ip1)/(2 * this->partial_x(i))) + ((theta * this->m_c(t_i, x_ip1))/pow(this->partial_x(i), 2));
        }            
    }

    return P;
}

Matrix PDE::computeQ(size_t i, size_t m, double theta)
{
    Matrix Q(m - 1, m - 1);
    double t_i = this->m_timeGrid[i];

    for(size_t j = 0; j < m - 1; j++)
    {
        double x_i = this->m_spaceGrid[j];
        double x_ip1 = this->m_spaceGrid[j+1];
        Q[j][j] = 1/this->partial_t(i) - (2 * (1 - theta) * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2);
        
        if(j + 1 < m - 1)
        {
            Q[j][j + 1] = ((1 - theta) * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2);
            Q[j + 1][j] = ((1 - theta) * this->m_c(t_i, x_ip1)) / pow(this->partial_x(i), 2);
        }            
    }

    return Q;
}

vector<double> PDE::computeV(size_t i, size_t m, double theta)
{
    vector<double> V(m - 1, 0.0);
    double t_i = this->m_timeGrid[i];

    for(int j = 0; j < m - 1; j++)
    {
        double x_i = this->m_spaceGrid[j];
        double x_ip1 = this->m_spaceGrid[j+1];
        if(j == 0)
            V[j] = this->m_d(t_i, x_i) + (-this->m_b(t_i, x_i)/(2*this->partial_x(i)) + (theta * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2)) * this->m_leftBoundary[i] + (((1 - theta) * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2)) * this->m_leftBoundary[i + 1];
        else if(j == m - 2)
            V[j] = this->m_d(t_i, x_i) + (-this->m_b(t_i, x_i)/(2*this->partial_x(i)) + (theta * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2)) * this->m_rightBoundary[i] + (((1 - theta) * this->m_c(t_i, x_i)) / pow(this->partial_x(i), 2)) * this->m_rightBoundary[i + 1];
        else
            V[j] = this->m_d(t_i, x_i);
    }

    return V;
}

// PDE Methods
void PDE::setPDE(function<double (double, double)>& a, function<double (double, double)>& b, function<double (double, double)>& c, function<double (double, double)>& d)
{
    this->m_a = a;
    this->m_b = b;
    this->m_c = c;
    this->m_d = d;
}

vector<double> PDE::createTerminalCondition(Option* option)
{
    vector<double> res;

    for(size_t i = 0; i < this->m_spaceGrid.size() - 1; i++)
    {
        res.push_back(option->getPayoff(this->m_spaceGrid[i]));
    }

    return res;
}

void PDE::setProblem(vector<double> leftBoundary, vector<double> rightBoundary, vector<double> terminalCondition)
{
    this->m_leftBoundary = leftBoundary;
    this->m_rightBoundary = rightBoundary;
    this->m_terminalCondition = terminalCondition;
}

void PDE::resolve()
{
    size_t m = this->m_spaceGrid.size();
    size_t n = this-> m_timeGrid.size();
    double theta = 0.5;
    
    //Set terminal condition
    this->m_Solution[n - 1] = this->m_terminalCondition;

    // Numerical propagation to solve the PDE
    for(int i = n - 2; i >= 0; i--)
    {
        Matrix P = computeP(i, m, theta);
        Matrix Q = computeQ(i, m, theta);
        vector<double> V = computeV(i, m, theta);
        vector<double> rhs = Q * this->m_Solution[i + 1] + V;
        this->m_Solution[i] = -P.invert() * rhs;
    }
}