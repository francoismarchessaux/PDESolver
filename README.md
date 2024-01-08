# M2 203 C++ Project
### _Group : Louis Faverjon, Camille Vallée, François Marchessaux_

## Project Goal and code implementation
The goal of this project is to implement a PDE solver to price a European Vanilla option.

To this end we have implemented several classes : 
- Class **_Option_**: This class allows to create a European vanilla option with its parameters, and to compute the payoff of the option given a certain spot with the function ```getPayoff()```.
- Class **_BlackScholes_**: This class uses the closed formula of the Black-Scholes model to compute the price of an option with the function ```price()```.
- Class **_Matrix_**: This class creates a Matrix object, and overrides mathematical operators to perform computations with matrices, vectors and scalars. The function ```invert()``` checks if a matrix is invertible, and if yes, inverts it with the Gaussian elimination method.
- Class **_PDE_**: This class allows to solve a PDE with the finite difference method. The constructor initalizes the time grid and space grid for the discretization of the problem, the terminal condition which is given by the payoff of the option at maturity, and the boundaries. The price of the option is bounded by 0 and by the maximum discounted spot of the space grid. Then, the function ```resolve()``` will use the finite difference method to solve the PDE. Finally, the price is extracted by the function ```solution()```, which performs a linear interpolation to find the price of the option at t = 0 and corresponding to the current spot price.

## Results

We have used these parameters for the option : 
| Parameter | Value |
| ------ | ------ |
|Option Type|Call|
|S|$100|
|K|$110|
|T|1Y|
|r|2%|
|σ|20%|
|div|0%|
|repo|0%|

By identification, we set the coefficient functions as follows : 
```math
a(t, x) = -r
```
```math
b(t, x) = x * (r - div - repo)
```
```math
c(t, x) = \frac{σ^2}{2} * x^2
```
```math
d(t, x) = 0
```

For the 1st simulation, we used these parameters for the PDE and obtained the following prices : 
| Parameter | Value |
| ------ | ------ |
|time steps|50|
|space steps|75|
|multiplier|4|
|**BS Price**|$4.944|
|**PDE Price**|$5.072|
|**Absolute Difference**|$0.128|

For the 2nd simulation, we decreased the time steps and obtained the following prices : 
| Parameter | Value |
| ------ | ------ |
|time steps|10|
|space steps|75|
|multiplier|4|
|**BS Price**|$4.944|
|**PDE Price**|$4.132|
|**Absolute Difference**|$0.812|

For the 3rd simulation, we decreased the space steps and obtained the following prices : 
| Parameter | Value |
| ------ | ------ |
|time steps|50|
|space steps|30|
|multiplier|4|
|**BS Price**|$4.944|
|**PDE Price**|$4.213|
|**Absolute Difference**|$0.731|

We can then observe a convergence of the price as time or space steps increase. 
