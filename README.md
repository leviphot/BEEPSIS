# BEEPSIS: Bayes estimation of experimental parameters in stochastic inertial systems

Functions in BEPPESIS package (beeepsis.m, beepsis_bin.m, and beepsis_ndbin.m) estimate the parameters influencing motion of a particle(s) in inertial system with random force. One can estimate the force profile, damping coefficient as well as the temperature of surrounding medium using a position only record of the particle trajectory. 

For details and examples of usage see https://arxiv.org/abs/2212.14043

Please note that the resulting "force" is normalized by particle mass while the "effective temperature" is given relative to $k_B T_{\mathrm{ref}}/m$ (where $T_{\mathrm{ref}}$ is some reference temperature and $m$ is particle mass).

## Installation

Download the BEEPSIS package and compile the C++ code using MATLAB *mex* compiler

```
        mex -R2018a solveBeepsisArg.cpp
```
Optionally, add BEEPSIS folder to MATLAB path

## Data preparation

Prior employing BEEPSIS prepare following data (variables)

1. **x**: particle trajectory (vector or matrix). Each column represents one dimension, rows are positions in time steps. Should be calibrated in meters.

2. **dt**: time step, 1 number

3. **fun**: function describing force acting on a particle in position *x* parametrized by a set of parameters *par*.\
 **Required only for parametric estimation using beepsis.m** \
One can use anonymous function as well. If force is time dependent, add extra "time" column to **x** variable. \
Function prototype is 
```MATLAB
        function F = fun(par, x)
```
4. **init**: Initial guesses of force parameters. Number / row vector / matrix\
One can test multiple guesses of the parameters. In such case each row of matrix is one initial guess of force parametrization\
**Required only for parametric estimation using beepsis.m** \
**Minimization procedure looks for scaling factors of the initial guess. Do not set the value to 0 (zero), algorithm would bot be able to obtain any other value**

5. **Gamma**: initial guess (or value) of damping parameter $\Gamma$\
can be column vector for beepsis.m in order to test multiple initial guesses

6. **kT_m**: ratio of temperature over particle mass multiplied by the Boltzmann constant (either initial guess, or value)\
**Required only for parametric estimation using beepsis.m** \
can be column vector for beepsis.m in order to test multiple initial guesses

7. **sigma**: value of standard deviation of the additive detection error of positions (optional). For multidimensional problems can be either one value or row vector with different values for each dimension


### Some examples of the parametrized force profiles
```MATLAB
        harmonic = @(omega, x) - mass*omega^2 * x; % harmonic oscillator
        duffing  = @(par, x)   - mass*par(1)^2 * x .* (1-par(2)*x.^2*1e12); % Duffing oscillator, Duffing strength is assumed to be in \mu m^{-2}, therefore 1e12

        %% force profile described by a spline with K knots
        K = 11; % number of knots
        knots_x      = linspace(min(x), max(x), K);  % knots positions
        knots_x(1)   = knots_x(1) - eps(knots_x(1)); % slightly shift the outermost knots
        knots_x(end) = knots_x(end) + eps(knots_x(end));
        spline_force = @(p, x) spline(knots_x, p, x); 
```

## Parametric force, damping and temperature estimation (**beepsis.m**)
For the estimation of force profile parameters, damping $\Gamma$ and effective temperature $T$ one should use function beepsis.m
```MATLAB
        function res = beepsis(x, dt, fun, init, Gamma, kBT_m, methods, sigma)
```
The meaning of all arguments (except *methods*) is described above. \
Detection error *sigma* is optional


**method** is a MATLAB char array (or string) that describes a sequence of minimization steps. In each step a different parameters may by estimated using different estimators and different minimization algorithm. Results of each step then may be used as input in a next sequence step. Or the different sequence may be started with the input values of parameters. 

**method** string has following structure

1. optional numbers *0-3* - for selection of algorithm
2. optional character *^* - to switch on parallel minimization
3. one of letters *vpcCGTA* for selection of parameters that are estimated and for selection of estimator 
4. steps 1 and 2 may be repeated to define a sequence that uses outputs of a previous step as an input in next one
5. optional character *,|* indicating start of a next sequence with initial guesses taken from the inputs of the **beepsis** function
6. new sequence definition from steps 1-3

An example of **method** string 
```MATLAB
        method = '0A|3^v1cg2^A|2^A';
```
The meaning of this example will be explained bellow. 

### Selection minimization function
character 0-3, optionally character '^' to use parallel computing capabilities of some of the minimization functions

1. **0**: no minimization is performed, only the value of estimator is returned. This is useful together with "multiple initial guesses" to map the landscape of the  the estimator
2. **1**: MATLAB function ***fminesearch*** is used for minimization. No extra options are applied
3. **2**: MATLAB function ***fminunc*** is used for minimization. (default value)
4. **2^**: MATLAB function ***fminunc*** is used for minimization. Parallel evaluation is switched on
5. **3**: MATLAB function **ga** - genetic algorithm (requires Global Optimization Toolbox). \
&nbsp;&nbsp;&nbsp;&nbsp;Lower and upper bound of force parameters is $\langle -4, 5\rangle \times$ the input value.\
&nbsp;&nbsp;&nbsp;&nbsp;Lower and upper bound of damping $\Gamma$ and temperature $T$ is $\langle 0.1, 10\rangle \times$ the input value.
6. **3^**: MATLAB function **ga** as above.  Parallel evaluation is switched on


### Selection of estimator and estimated parameters
One of letters 'vpcCGTA' which define an estimator formula and a set of parameters which are used in minimization 

1. **v**: velocity only dependent PDF (Eq. B1); minimization for force; damping and temperature are fixed
2. **p**: position and velocity dependent PDF (Eq. B2); minimization for force; damping and temperature are fixed
3. **c**: position only dependent PDF (Eqs. 7-11); minimization for force; damping and temperature are fixed  
4. **C**: position only dependent PDF (Eqs. 7-11); minimization for force and damping; temperature is fixed  
5. **G**: position only dependent PDF (Eqs. 7-11); minimization for damping; force and temperature are fixed  
6. **T**: position only dependent PDF (Eqs. 7-11); minimization for damping and temperature;  force is fixed  
7. **A**: position only dependent PDF (Eqs. 7-11); minimization for force, damping, and temperature


### Multi step minimization
It is often desired to perform parameter search in several steps with results of one step taken as inputs to another. This approach usually leads to more "stable" solutions and avoids un-physical values of searched parameters.  In the case of BEEPSIS pone may want to minimize for force parameters first, then damping and temperature and finally for all parameters close to already found minimum of estimator. 

In case of BEEPSIS in is possible to chain several definitions of minimization functions and  estimators/parameters with results of one step being taken as inputs of a next step. 

One may also want to try several possible chains with each starting in the initial guess. For such separation of chains use one of characters **,|**

Lets now look back to the previous example
```MATLAB
        method = '0A|3^v1cg2^A|2^A';
```
It defines 3 separate chains 

1. **0A** no minimization, just check value of estimator with position only dependent PDF. 
   
2. 
    1. **3^v** use genetic algorithm in parallel to look for parameters of force using velocity only dependent PDF (Eq. B1);
   
    2. **1c** use ***fminsearch*** to look for parameters of force
    3. **g** use ***fminunc*** to look for damping
    4. **2^A** use ***fminunc*** to look for all force, damping, and temperature. Use parallel computations
   
3. **2^A** use ***fminunc*** to look for all force, damping, and temperature. Start minimization directly from the initial guess. Use parallel computations


### Multiple initial guesses
In order to test several initial guesses of the parameters (force, damping, temperature) **beepsis.m** accepts vectors (matrix) of these parameters. However some restrictions apply:

1. force function parameters: **each row is a new initial guess, multiple parameters must be given in columns**
2. multiple initial guesses of damping and temperature could be either row or column vectors
3. **number of initial guesses for force, damping and temperature must be either the same or one**. If only guess is supplied (e.g. for temperature) this value is used in combination with other initial guesses

The **method** chain is then applied to all of the initial guesses


### Results
Results of **beepsis.m** are returned in a structure (or array of structures). Fields of this structure are vectors (arrays) with values obtained in each minimization step. 

fields of **res** structure are

| | | 
|---|---|
|**res.force**|  column vector of array with each row corresponding to the value of force function parameters obtained  by the minimization |
|**res.Gamma**|  column vector of minimized damping coefficients|
|**res.effTr**|  column vector of effective relative temperatures, i.e. values are factors which multiply input **kBT_m**|
|**res.forceerr** <br>  **res.Gammaerr** <br>**res.effTrerr**| estimates of error of the parameters based on inversion of Hessian matrix.  Values are given with 95% confidence interval. <br> If the minimization for a given value is not performed, then error is NaN  and the returned value is either input or value from the previous step |
|**res.covmat**| cell array, covariance matrix for the minimized  parameters in given step |
|**res.Svalue**| value of the minimized $-\log(P)$|
|**res.method**| minimization method of current step (letters from the above definition)
|**res.algorithm**|  function used for minimization, one of numbers 0-3 (see above)|


## Nonparametric force estimation in 1D (**beepsis_bin.m**)
For the estimation of unknown force profile parameters in 1D one should use function **beepsis_bin.m**
```MATLAB
         function res = beepsis_bin(x, dt, bins,  Gamma, methods)
```
In contrast to **beepsis.m** no minimization is performed, but one can estimate only force profiles, and it is assumed that the detection error is negligible. 
The theoretical background is described in Appendix E of the main text

### Inputs
|||
|---|---|
|**x**| particle trajectory vector, should be in meters|
|**dt**| time step, 1 number | 
|**bins**|  number of bins/bin edges<br> empty     bin count automatically calculated <br>  1 number  maximal number of bins (could be adjusted, see bellow) <br> vector of bin edges  |
|**Gamma**| Damping coefficient $\Gamma$ in [s$^{-1}$] |
|**methods**| simplified character string defining method of force  estimation. <br> Can be any combination of following<br> 'v'  inference for velocities,according to Eq. (E9) <br> 'p'  inference for positions and init velocity according to Eq. (E10) <br>'c'  inference for time correlated misfits in positions,            according to Eqs (E4-E8)|

### Internal settings variables

|  |  | 
|---|---|
|**NMIN**     = 20| minimal count of data points in one bin. In the count is lower adjacent bins are merged|
|**MAXBIN**    = 100| maximal number of bins|
|**CINVTERMS** = 25| number of "upper and lower diagonals" taken into calculation using covariance matrix inversion by Eqs. (E7-E8)|

### Output
Results are filed in **res** structure as follows

| | |
|---|---|
|**res.bincenter**|centers of each bin (useful for plotting)|
|**res.binedges**| edges of bins|
|**res.forcev**| force calculated using Eq. (E9), if **method** string contained **v**|
|**res.forcep**| force calculated using Eq. (E9), if **method** string contained **p**|
|**res.force**| force calculated using Eq. (E4-E8), if **method** string contained **c**|
|**res.vdrift**| drift velocity averaged in each bin |

## Nonparametric force estimation in multiple dimensions (**beepsis_ndbin.m**)
For the estimation of unknown force profile parameters in  multiple dimensions. Extension of **beepsis_bin.m**
```MATLAB
         function res = beepsis_ndbin(data, dt, nbin, Gamma, methods)
```

### Inputs

| | | 
|---|---|
|**data**| matrix of positions (rows are time, columns degrees of freedom), in meters|
|**dt**| time step, 1 number | 
|**bins**|  number of bins<br> empty: bin count automatically calculated <br>  1 number:  number of bins, same for all degrees of freedom <br> bin counts for, each number for 1 degree of freedom |
|**Gamma**| Damping coefficient $\Gamma$ in [s$^{-1}$] |
|**methods**| simplified character string defining method of force  estimation. <br> Can be any combination of following<br> 'v'  inference for velocities,according to Eq. (E9) <br> 'p'  inference for positions and init velocity according to Eq. (E10) <br>'c'  inference for time correlated misfits in positions,            according to Eqs (E4-E8)|

### Internal settings variables

||  |
|---|---|
|**MAXBIN**    = 100| maximal number of bins|
|**CINVTERMS** = 25| number of "upper and lower diagonals" taken into calculation using covariance matrix inversion by Eqs. (E7-E8)|

### Output
Results are filed in **res** structure as follows

| | |
|---|---|
|**res.bincenter**|centers of each bin (useful for plotting)|
|**res.binedges**| edges of bins|
|**res.forcev**| force calculated using Eq. (E9), if **method** string contained **v**|
|**res.forcep**| force calculated using Eq. (E9), if **method** string contained **p**|
|**res.force**| force calculated using Eq. (E4-E8), if **method** string contained **c**|
|**res.vdrift**| drift velocity averaged in each bin |




