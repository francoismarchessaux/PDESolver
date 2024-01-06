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

    // PDE Parameters
    size_t timeSteps = 50;
    size_t spaceSteps = 100;
    double multiplier = 10;
    function<double(double, double)> a = [r](double t, double x) { return -r; };
    function<double(double, double)> b = [r](double t, double x) { return r * x; };
    function<double(double, double)> c = [vol](double t, double x) { return (pow(vol, 2) / 2) * pow(x, 2); };
    function<double(double, double)> d = [](double t, double x) { return 0; };

    // Create Vanilla European option
    Option option(T, K, OptionType::Put);
    cout << option << endl;
    
    // Compute option price with BlackScholes
    BlackScholes bs(option, S, r, vol);
    cout << "Price of the option with Black-Scholes is $" << bs.price() << endl;

    // Compute option price with PDE
    PDE pde(option, timeSteps, T, spaceSteps, multiplier, vol, S, r, a, b, c, d);
    pde.resolve();
    cout << "Price of the option with PDE is $" << pde.solution(S) << endl;

    return 0;
}
