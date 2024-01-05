#include <iostream>
#include <functional>
#include "Matrix.h"
#include "Option.h"
#include "BlackScholes.h"
#include "PDE.h"

using namespace std;

int main()
{
    // Options Parameters
    double S = 100.0;
    double K = 110.0;
    double T = 1.0;
    double r = 0.02;
    double vol = 0.2;

    // Create Vanilla European Call option
    Option callOption(T, K, OptionType::Put);
    cout << callOption << endl;
    
    // Compute option price with BlackScholes
    BlackScholes bs_pricer(&callOption, S, r, vol);
    cout << "Price of the option with Black-Scholes is $" << bs_pricer() << endl;

    // PDE parameters
    size_t timeSteps = 50;
    double minTime = 0.0;
    double maxTime = T;
    size_t spaceSteps = 100;
    double multiplier = 10;
    PDE pde(timeSteps, minTime, maxTime, spaceSteps, multiplier, vol, S);
    vector<double> spaceGrid = pde.getSpaceGrid();
    vector<double> timeGrid = pde.getTimeGrid();
    vector<double> leftBoundary(timeSteps, 0.0);
    double maxSpot = *max_element(spaceGrid.begin(), spaceGrid.end());
    vector<double> rightBoundary;

    for(int i = 0; i < timeGrid.size(); i++)
    {
        rightBoundary.push_back(exp(-r * (timeGrid[timeGrid.size()] - timeGrid[i])) * maxSpot);
    }

    // Coefficient functions
    function<double(double, double)> a;
    a = [r](double t, double x) {
        return -r;
    };
    function<double(double, double)> b;
    b = [r](double t, double x) {
        return r*x;
    };
    function<double(double, double)> c;
    c = [vol](double t, double x) {
        return (pow(vol, 2) / 2) * pow(x, 2);
    };
    function<double(double, double)> d;
    d = [](double t, double x) {
        return 0;
    };

    // Set up PDE
    vector<double> terminalCondition = pde.createTerminalCondition(&callOption);
    pde.setProblem(leftBoundary, rightBoundary, terminalCondition);
    pde.setPDE(a, b, c, d);

    // Solve PDE
    pde.resolve();

    // Display PDE solution
    double pde_price = pde.solution(S);
    cout << "Price of the option with PDE is $" << pde_price << endl;

    return 0;
}
