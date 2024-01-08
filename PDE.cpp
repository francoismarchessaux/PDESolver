#include <cmath>
#include <algorithm>
#include "Option.h"
#include "Matrix.h"
#include "PDE.h"

// Overriding + operator to sum two vectors
std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2) 
{
    if (v1.size() != v2.size())
        throw std::invalid_argument("Vector addition: vectors must have the same size");

    std::vector<double> res(v1.size());

    for (size_t i = 0; i < v1.size(); ++i)
        res[i] = v1[i] + v2[i];

    return res;
}

// Initializer Constructor
PDE::PDE(Option option, size_t nTimeSteps, double T, size_t nSpaceSteps, double multiplier, double vol, double spot, double r, std::function<double (double, double)>& a, std::function<double (double, double)>& b, std::function<double (double, double)>& c, std::function<double (double, double)>& d) : m_Option(option), m_a(a), m_b(b), m_c(c), m_d(d)
{
    // Initialize the time grid
    setTimeGrid(nTimeSteps, T);

    // Initialize the space grid
    setSpaceGrid(nSpaceSteps, multiplier, vol, spot);

    // Initialize solution matrix
    m_Solution = Matrix(nTimeSteps - 1, nSpaceSteps - 1);

    // Set the terminal condition
    setTerminalCondition();

    // Set the left and right boundaries
    setBoundaries(nTimeSteps, r);
}

// Accessors
std::vector<double> PDE::getTimeGrid() const { return m_timeGrid; }
std::vector<double> PDE::getSpaceGrid() const { return m_spaceGrid; }

// Mutators
void PDE::setTimeGrid(size_t nTimeSteps, double T)
{
    double deltaTime = T / static_cast<double>(nTimeSteps);

    for(size_t i = 0; i < nTimeSteps; i++)
        m_timeGrid.push_back(i * deltaTime);
}

void PDE::setTimeGrid(const std::vector<double>& timeGrid)
{
    m_timeGrid = timeGrid;
    sort(m_timeGrid.begin(), m_timeGrid.end());
}

void PDE::setSpaceGrid(size_t nSpaceSteps, double multiplier, double vol, double spot)
{
    double lowerBoundary = (spot - multiplier * vol * spot) >= 0 ? spot - multiplier * vol * spot : 0;
    double upperBoundary = spot + multiplier * vol * spot;
    double deltaSpace = (upperBoundary - lowerBoundary) / static_cast<double>(nSpaceSteps);

    for(size_t i = 0; i < nSpaceSteps; i++)
        m_spaceGrid.push_back(lowerBoundary + i * deltaSpace);
}

void PDE::setSpaceGrid(const std::vector<double>& spaceGrid)
{
    m_spaceGrid = spaceGrid;
    sort(m_spaceGrid.begin(), m_spaceGrid.end());
}

void PDE::setTerminalCondition()
{
    // Terminal condition is the payoff at maturity
    std::vector<double> res;
    for(size_t i = 0; i < m_spaceGrid.size() - 1; i++)
        res.push_back(m_Option.getPayoff(m_spaceGrid[i]));

    m_terminalCondition = res;
}

void PDE::setBoundaries(size_t nTimeSteps, double r)
{
    // Minimum price is bounded by 0
    std::vector<double> leftBoundary(nTimeSteps, 0.0);

    // Maximum price is bounded by max spot actualized with remaining maturity
    std::vector<double> rightBoundary;
    double maxSpot = *max_element(m_spaceGrid.begin(), m_spaceGrid.end());
    double div = m_Option.getDiv();
    double repo = m_Option.getRepo();
    for(int i = 0; i < nTimeSteps; i++)
        rightBoundary.push_back(exp((r - div - repo) * (m_timeGrid[nTimeSteps - 1] - m_timeGrid[i])) * maxSpot);

    m_leftBoundary = leftBoundary;
    m_rightBoundary = rightBoundary;
}

// Partial derivatives of t and x
double PDE::partial_t() const { return m_timeGrid[1] - m_timeGrid[0]; }
double PDE::partial_x() const { return m_spaceGrid[1] - m_spaceGrid[0]; }

// Matrices computations
Matrix PDE::computeP(size_t i, size_t m, double theta)
{
    Matrix P(m - 1, m - 1);
    double t_i = m_timeGrid[i];

    for(size_t j = 0; j < m - 1; j++)
    {
        double x_i = m_spaceGrid[j];
        double x_ip1 = m_spaceGrid[j+1];
        P[j][j] = m_a(t_i, x_i) - ((1.0 / partial_t()) + ((2.0 * theta * m_c(t_i, x_i)) / pow(partial_x(), 2)));
        
        if(j + 1 < m - 1)
        {
            P[j][j + 1] = (m_b(t_i, x_i)/(2 * partial_x())) + ((theta * m_c(t_i, x_i))/pow(partial_x(), 2));
            P[j + 1][j] = -(m_b(t_i, x_ip1)/(2 * partial_x())) + ((theta * m_c(t_i, x_ip1))/pow(partial_x(), 2));
        }            
    }

    return P;
}

Matrix PDE::computeQ(size_t i, size_t m, double theta)
{
    Matrix Q(m - 1, m - 1);
    double t_i = m_timeGrid[i];

    for(size_t j = 0; j < m - 1; j++)
    {
        double x_i = m_spaceGrid[j];
        double x_ip1 = m_spaceGrid[j+1];
        Q[j][j] = 1/partial_t() - (2 * (1 - theta) * m_c(t_i, x_i)) / pow(partial_x(), 2);
        
        if(j + 1 < m - 1)
        {
            Q[j][j + 1] = ((1 - theta) * m_c(t_i, x_i)) / pow(partial_x(), 2);
            Q[j + 1][j] = ((1 - theta) * m_c(t_i, x_ip1)) / pow(partial_x(), 2);
        }            
    }

    return Q;
}

std::vector<double> PDE::computeV(size_t i, size_t m, double theta)
{
    std::vector<double> V(m - 1, 0.0);
    double t_i = m_timeGrid[i];

    for(int j = 0; j < m - 1; j++)
    {
        double x_i = m_spaceGrid[j];
        double x_ip1 = m_spaceGrid[j+1];
        if(j == 0)
            V[j] = m_d(t_i, x_i) + (-m_b(t_i, x_i)/(2*partial_x()) + (theta * m_c(t_i, x_i)) / pow(partial_x(), 2)) * m_leftBoundary[i] + (((1 - theta) * m_c(t_i, x_i)) / pow(partial_x(), 2)) * m_leftBoundary[i + 1];
        else if(j == m - 2)
            V[j] = m_d(t_i, x_i) + (-m_b(t_i, x_i)/(2*partial_x()) + (theta * m_c(t_i, x_i)) / pow(partial_x(), 2)) * m_rightBoundary[i] + (((1 - theta) * m_c(t_i, x_i)) / pow(partial_x(), 2)) * m_rightBoundary[i + 1];
        else
            V[j] = m_d(t_i, x_i);
    }

    return V;
}

// PDE Methods
void PDE::resolve()
{
    size_t m = m_spaceGrid.size();
    size_t n = m_Solution.getNRows();
    double theta = 0.5;
    
    //Set terminal condition
    m_Solution[n - 1] = m_terminalCondition;

    // Numerical propagation to solve the PDE
    for(int i = n - 2; i >= 0; i--)
    {
        Matrix P = computeP(i, m, theta);
        Matrix Q = computeQ(i, m, theta);
        std::vector<double> V = computeV(i, m, theta);
        std::vector<double> rhs = Q * m_Solution[i + 1] + V;
        m_Solution[i] = -P.invert() * rhs;
    }
}

double PDE::solution(double spot)
{
    // Get closest indexes of spot in space grid
    size_t iLower = 0;
    while(m_spaceGrid[iLower] <= spot)
        iLower++;

    // Linear interpolation to find the price
    double interpolated_price = m_Solution[0][iLower] + (spot - m_spaceGrid[iLower]) * ((m_Solution[0][iLower + 1] - m_Solution[0][iLower])/(m_spaceGrid[iLower + 1] - m_spaceGrid[iLower]));

    return interpolated_price;
}