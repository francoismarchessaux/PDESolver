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
    double PDEvol = 0.5;
    vector<double> leftBoundary(timeSteps, 0.0);
    vector<double> rightBoundary(timeSteps, S);

    function<double(double, double)> a;
    a = [r](double x, double y) {
        return -r;
    };
    function<double(double, double)> b;
    b = [r](double x, double y) {
        return r*x;
    };
    function<double(double, double)> c;
    c = [vol](double x, double y) {
        return (pow(vol, 2) / 2) * pow(x, 2);
    };
    function<double(double, double)> d;
    d = [](double x, double y) {
        return 0;
    };

    // Set up PDE
    PDE pde(timeSteps, minTime, maxTime, spaceSteps, multiplier, PDEvol, K);
    vector<double> terminalCondition = pde.createTerminalCondition(&callOption);
    pde.setProblem(leftBoundary, rightBoundary, terminalCondition);
    pde.setPDE(a, b, c, d);

    // Solve PDE
    pde.resolve();

    // Display PDE solution
    Matrix solutionMatrix = pde.getSolution();
    for(int i = 0; i < solutionMatrix.getNRows(); i++)
    {
        for(int j = 0; j < solutionMatrix.getNCols(); j++)
        {
            cout << solutionMatrix[i][j] << " ";
        }
        cout << endl;
    }

    // SOME TESTS
    /*
    vector<double> spaceGrid = pde.getSpaceGrid();
    vector<double> timeGrid = pde.getTimeGrid();
    
    cout << "\nTerminal Condition : " << endl;
    for(int i = 0; i < terminalCondition.size(); i++)
    {
        cout << terminalCondition[i] << " ";
    }
    
    cout << "\nTime Grid : " << endl;
    for(int i = 0; i < timeGrid.size(); i++)
    {
        cout << timeGrid[i] << " ";
    }

    cout << "\nSpace Grid : " << endl;
    for(int i = 0; i < spaceGrid.size(); i++)
    {
        cout << spaceGrid[i] << " ";
    }
    */
    //vector<vector<double> > values = {{1, 0, 2, -1, 3}, {3, 0, 0, 5, 2}, {2, 1, 4, -3, 1}, {1, 0, 5, 0, 4}, {3, 1, 6, 5, 2}};
    //Matrix mat(values);
    //mat.display();
    //cout << endl;

    return 0;
}
