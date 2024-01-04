#include <iostream>
#include <functional>
#include "Matrix.h"
#include "Option.h"
#include "BlackScholes.h"
#include "PDE.h"

using namespace std;

#if __cplusplus < 201103L
#error "C++11 or later is required."
#endif

double generateRandomDouble(double min, double max) 
{
    return min + static_cast<double>(rand()) / (RAND_MAX / (max - min));
}

// Function to create a vector<vector<double>> of size n filled with random values
std::vector<std::vector<double>> createRandomMatrix(int n, double minValue, double maxValue) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    // Seed for random number generation
    std::srand(static_cast<unsigned>(std::time(0)));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = generateRandomDouble(minValue, maxValue);
        }
    }

    return matrix;
}

int main()
{
    // Options Parameters
    double S = 100.0;
    double K = 110.0;
    double T = 1.0;
    double r = 0.02;
    double vol = 0.2;

    // Create Vanilla European option
    Option callOption(T, K, OptionType::Call);
    cout << callOption << endl;
    
    // Compute option price with BlackScholes
    BlackScholes bs_pricer(&callOption, S, r, vol);
    cout << "Price of the option with Black-Scholes is $" << bs_pricer() << endl;

    // PDE parameters
    size_t timeSteps = 50;
    double minTime = 0.0;
    double maxTime = 1.0;
    size_t spaceSteps = 100;
    double multiplier = 2;
    vector<double> leftBoundary(timeSteps, 0.0);
    vector<double> rightBoundary(timeSteps, S);

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
    PDE pde(timeSteps, minTime, maxTime, spaceSteps, multiplier, vol, K);
    vector<double> terminalCondition = pde.createTerminalCondition(&callOption);
    pde.setProblem(leftBoundary, rightBoundary, terminalCondition);
    pde.setPDE(a, b, c, d);

    // Solve PDE
    pde.resolve();

    // Display PDE solution
    Matrix solutionMatrix = pde.getSolution();
    vector<double> solution_t0 = solutionMatrix[0];
    for(size_t i = 0; i < solution_t0.size(); i++)
    {
        cout << solution_t0[i] << " ";
    }

    return 0;
}
