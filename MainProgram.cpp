#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include "Option.h"
#include "BlackScholes.h"
#include "Matrix.h"
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
    double div = 0.0;
    double repo = 0.0;

    // PDE Parameters
    size_t timeSteps = 50;
    size_t spaceSteps = 25;
    double multiplier = 4;

    // Coefficient functions found by identification
    function<double(double, double)> a = [r](double t, double x) { return -r; };
    function<double(double, double)> b = [r, div, repo](double t, double x) { return (r - div - repo) * x; };
    function<double(double, double)> c = [vol](double t, double x) { return (pow(vol, 2) / 2) * pow(x, 2); };
    function<double(double, double)> d = [](double t, double x) { return 0; };

    // Create Vanilla European option
    Option option(T, K, div, repo, OptionType::Call);
    cout << option << endl;
    
    // Compute option price with BlackScholes
    BlackScholes bs(option, S, r, vol);
    cout << "Price of the option with Black-Scholes is $" << std::fixed << std::setprecision(3) << bs.price() << endl;

    // Compute option price with PDE
    PDE pde(option, timeSteps, T, spaceSteps, multiplier, vol, S, r, a, b, c, d);
    pde.resolve();
    cout << "Price of the option with PDE is $" << std::fixed << std::setprecision(3) << pde.solution(S) << endl;

    return 0;
}

// BS Price : $4.944
// T steps : 50
// S steps : 75
// multiplier : 4
// Price : $5.072

// T steps : 10
// S steps : 75
// multiplier : 4
// Price : $4.132

// T steps : 50
// S steps : 30
// multiplier : 4
// Price : $4.213