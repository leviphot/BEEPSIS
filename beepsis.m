function res = beepsis(data, dt, fun, init, Gamma, kBT_m, methods, sigmax)
% function res = beepsis(data, dt, fun, init, Gamma, kBT_m, methods, sigmax)
% 
% evaluates normalized force profile, and/or damping, and/or medium effective temperature
% for the particle with inertia from its position  using Bayesian inference approach
% 
% INPUT
%   data     vector of positions (1D) or matrix (each column is one axis), should be in meters
%   dt       time step between measured positions
%   fun      function that calulates normalized force (force/mass) assuming 
%            parcile is in x0
%            function must accept arguments in a shape fun(par, x)
%               where par is a vector of parameters (same size as init)
%               and x is a vector (!D case) or matrix for multidimensional case
%   init     initialization guess of the normalized force profile parametrization
%            should be row vector
%            if init is colum or matrix then minimization is performed for each row as
%            a staring point (i.e. multiple starting points)
%   Gamma    damping coefficient or its initial guess
%            multiple starting points minimization is performed if Gamma is column vector
%            combines with rows in init and values in kT_m
%   kT_m     Boltzman constant times environment temperature over mass
%            multiple starting points minimization is performed for each value as a
%            starting point, combines with rows in init and values in Gamma
%   methods  character string defing method of force (and other) parameter
%            estimation. It can define several chain of several consequent
%            minimizations for different variables using different
%            minmizing functions. 
%            Can be any combination of following
%            '0-3'  for selection of minimizer method (function)
%               '0'  no minimization, just calculates value of -log(likehood) for given
%                    parameters
%               '1'  use fminsearch
%               '2'  use fminunc (default)
%               '3'  use genetic algorithm ga, ranges for parameterer search
%                  are set to +-5 times initial guess
%             letters 'vpcVPCGTA' defining which an approach and parameter
%               for minimization
%               'v'  velocitity only dependent PDF (Eq. B1)
%                    minimization for force, damping and temperature fixed
%               'p'  positon and velocty dependent PDF (Eq. B2)
%                    minimization for force, damping and temperature fixed
%               'c'  position only dependent PDF (Eqs. 7-11)
%                    minimization for force, damping and temperature fixed 
%               'C'  same as small 'c'
%                    minimization for force and damping, T fixed
%               'G'  same as small 'c'
%                    minimization for damping only, force and T fixed
%               'T'  same as small 'c'
%                    minimization for damping and temperaure, force fixed
%               'A'  same as 'c' 
%                    minimization for force, damping and T
%             the combinantion of optional number and one of above
%               mentioned letters define a chain of consequent
%               minimizations with output of previous entering as an input
%               to the next one
%             space or ',' defines a new chain of minimizations starting
%               again with the input parameters
%   sigmax   estimate of white noise sigma of position uncertainty in meters
%            can be vector for multidimensional input data
%
% OUTPUT
%   res      structure (or structure array) with results of minimizations 
%            each field is a vector (or array) representing one step in
%            chain minimization so one can see the progress
%            
%            res.force  column vector of array with each row corresponding
%                       to the value of force function parameters obtained 
%                       by the minimization 
%            res.Gamma  column vector of minimized damping coefficients
%            res.effTr  column vector of effective relative temperatures, 
%                       i.e. values are factors which multiply input kBT_m
%            res.forceerr
%            res.GammaErr
%            res.effTr  estimates of error of the parameters
%                       if the value is not minimized, then error is NaN
%                       and the returned value is eother input or value
%                       from the previous step
%            res.covmat cell array, covariance matrix for the mimized
%                       parameters in given step
%            res.Svalue value of the minimized -log(P)
%            res.method minimizarion method of current step (letters from
%                       above definition)
%            res.algorithm  function used for minimization, one of numbers
%                       0-3 (see above)
%
% This file is part of the BEEPSIS toolbox.
% See LICENSE.md for information about using/distributing this file.


if nargin<=6
    methods = 'C';
end


if nargin<=7
    sigmax = 0; % standard deviation of detection error
end

% multiple inputs 
if numel(Gamma)>1 || numel(kBT_m) > 1 || size(init,1) > 1
    NG = numel(Gamma);
    NT = numel(kBT_m);
    NI = size(init,1);
    pc = [NI, NG, NT];
    pc = pc(pc>1);
    if numel(pc)> 1 || any(pc(2:end)-pc(1) ~= 0)
        error('BEEPSIS:multipleInits:inconsistentSize', 'BEEPSIS inconsitent sizes of initial guesses for multiple initial guesses')
    end
    for kk = 1:pc(1)
        if NI == 1
            il = init;
        else
            il = init(kk,:);
        end
        if NG == 1
            ig = Gamma;
        else
            ig = Gamma(kk);
        end
        if NT == 1
            ik = kBT_m;
        else
            ik = kBT_m(kk);
        end
        rl = beepsis(data, dt, fun, il, ig, ik, methods, sigmax);
        if kk == 1
            res = rl;
            res = repmat(res, pc(1), 1);
        else
            res(kk) = rl;            
        end
    end
    return
end

chain = analyze_beepsis_method(methods, numel(init));

x0  = data(2:end-1,:);
x1  =  data(3:end,:);
xm1 =  data(1:end-2,:);

NX  =  numel(x0);


GT    = Gamma*dt;
EGT   = exp(-GT);
M1EGT = -expm1(-GT);

res.force     = NaN(numel(chain), numel(init));
res.forceerr  = NaN(numel(chain), numel(init));
res.Gamma     = NaN(numel(chain), 1);
res.Gammaerr  = NaN(numel(chain), 1);
res.effTr     = NaN(numel(chain), 1);
res.effTrerr  = NaN(numel(chain), 1);
res.covmat    = cell(numel(chain), 1);
res.Svalue    = NaN(numel(chain), 1);
res.method    = char(numel(chain), 1);
res.algorithm = NaN(numel(chain), 1);

for kk = 1:numel(chain)
    
    if kk > 1 && chain(kk).useprevsol
        finit = res.force(kk-1,:);
        ginit = res.Gamma(kk-1);
        tinit = res.effTr(kk-1);
    else
        finit = init;
        ginit = Gamma;
        tinit = 1;
    end
    abspar = 0;
    switch chain(kk).method
        case 'v'
            if all(sigmax == 0)
                sfun  = @(scale) sum(((x1-x0) - (x0-xm1)*EGT - fun(scale.*finit, x0)*dt/Gamma*M1EGT).^2, 'all');
            else
                sfun  = @(scale) gsumfun_v(scale.*finit, ginit, x1, x0, xm1, fun, dt, kBT_m, sigmax);
            end
            sinit = finit;
        case 'p'
            if all(sigmax == 0)
                sfun  = @(scale) sum((x1-x0 - (x1-xm1)/2/GT*M1EGT - fun(scale.*finit, x0)/Gamma^2*(GT - M1EGT)).^2, 'all');
            else
                sfun  = @(scale) gsumfun_p(scale.*finit, ginit, x1, x0, xm1, fun, dt, kBT_m, sigmax);
            end
            sinit = finit;
        case 'c'
            sfun  = @(scale) gsumfun_c(scale.*finit,  Gamma, x1, x0, xm1, fun, dt, kBT_m, sigmax);
            sinit = finit;
        case 'C'
            sfun  = @(scale) gsumfun_c(scale(1:end-1).*finit, scale(end)*ginit, x1, x0, xm1, fun, dt, kBT_m, sigmax);
            sinit = [finit, ginit];
        case 'G'
            sfun  = @(scale) gsumfun_c(finit, scale(1)*ginit, x1, x0, xm1, fun, dt, kBT_m, sigmax);
            sinit = ginit;
        case 'T'
            sfun  = @(scale) gsumfun_c(finit, scale(1)*ginit, x1, x0, xm1, fun, dt, scale(2)*kBT_m, sigmax);
            sinit = [ginit, tinit];
        case 'A'        
            sfun  = @(scale) gsumfun_c(scale(1:end-2).*finit, scale(end-1)*ginit, x1, x0, xm1, fun, dt, scale(end)*kBT_m, sigmax);
            sinit = [finit, ginit, tinit];
    end

    if chain(kk).minfun == 0
        fval = sfun(ones(size(sinit)));
        p  = sinit;
        e  = NaN(size(sinit));
        CM = NaN;
    else
        [p, e, CM, fval] = fit_likehood(sfun, sinit, NX, abspar, [chain(kk).minfun, chain(kk).parallel]);        
    end

    if all(chain(kk).method ~= 'GT')
        res.force(kk,:)     = p(1:numel(init));
        res.forceerr(kk,:)  = e(1:numel(init));
    else
        res.force(kk,:)     = finit;
    end
    if any(chain(kk).method == 'vpc')
        res.Gamma(kk) = ginit;
    elseif any(chain(kk).method == 'TA')
        res.Gamma(kk)    = p(end-1);        
        res.Gammaerr(kk) = e(end-1);
    else
        res.Gamma(kk)    = p(end);        
        res.Gammaerr(kk) = e(end);
    end
    if any(chain(kk).method == 'TA')
        res.effTr(kk)    = p(end);
        res.effTrerr(kk) = e(end);
    else
        res.effTr(kk) = tinit;
    end
    res.covmat{kk}    = CM;
    res.Svalue(kk)    = fval;
    res.method(kk)    = chain(kk).method;
    res.algorithm(kk) = chain(kk).minfun;

end
        
end


function chain = analyze_beepsis_method(method, parameter_count)
    if isempty(method)
        method = '1c';
    end
    if ~(ischar(method) || isstring(method))
        error('BEEPSIS:methodArray', 'BEEPSIS method is not char array')
    end
    valid_method = 'vpcCGTA';
    method_par   = [0 0 0 1 1 1 -1 -2 2];
    valid_fun    = '0123';
    valid_sep    = ',-|';
    valid_par    = '^';

    ichain = 1;
    ac.method     = '';
    ac.minfun     = 2;
    ac.useprevsol = false;
    ac.par_count  = parameter_count;
    ac.parallel   = false;
    chain = struct(ac);

    for kk = 1:numel(method)
        if any(method(kk) == valid_sep)
            ac.useprevsol = false;
            continue;
        end
        if any(method(kk) == valid_fun)
            ac.minfun = str2double(method(kk));
        end
        if any(method(kk) == valid_par)
            ac.parallel = true;
        end
        if any(method(kk) == valid_method)
            ac.method = method(kk);
            ii = find(method(kk) == valid_method);
            if method_par(ii) < 0
                ac.par_count = abs(method_par(ii));
            else
                ac.par_count = parameter_count + method_par(ii);
            end
            chain(ichain) = ac;
            ichain = ichain +1;
            ac.useprevsol = true;
        end
    end
end


function [par, err, CM, fval] = fit_likehood(sfun, init, NX, abs_from_endc, minmethod)
    par = ones(size(init));
    hessian = [];
    if numel(minmethod) == 2
        useparallel = minmethod(2);
        minmethod = minmethod(1);
    end
    if minmethod == 1
        par = fminsearch(sfun, par);
    elseif minmethod == 2
        options = optimoptions(@fminunc,'Display','off', 'UseParallel', logical(useparallel));
        [par, fval,~,~,~,hessian] = fminunc(sfun, par, options);
    elseif minmethod == 3
        ub = ones(1,numel(par)) * 5;
        lb = ones(1,numel(par)) * -4;
        if abs_from_endc == 1
            lb(end) = 0.1;
            ub(end) = 10;
        elseif abs_from_endc == 2
            lb(end-1) = 0.1;
            ub(end-1) = 10;
        end
        
        options = optimoptions(@ga,'Display','off', 'UseParallel', logical(useparallel));
        
        par = ga(sfun, numel(par), [],[],[],[],lb, ub, [], options);
        hessian = [];
    end
    if isempty(hessian)
        hessian = numhessian(sfun, par);
        fval    = sfun(par);
    end
    
    vv = NX-numel(par);
    CM = inv(hessian);
    
    dc   = sqrt(diag(CM));
    err  = (dc) .* tinv(1-0.95/2,vv);
    
    dc    = diag(1./dc);
    CM = dc*CM*dc; %#ok<MINV>
    
    par = par.*init;
    err = err'.*init;  
    if nargin > 3 && abs_from_endc > 1
        par(end-(abs_from_endc-1):end) = abs(par(end-(abs_from_endc-1):end));
        err(end-(abs_from_endc-1):end) = abs(err(end-(abs_from_endc-1):end));
    end

end


function r = gsumfun_v(p, GM, x1, x0, xm1, fun, dt, kBT_m, sigma)
    A1 = numel(x0)/2;
    B  = kBT_m*2;
    
    GM    = abs(GM);
    GT    = GM*dt;
    EGT   = exp(-GT);
    M1EGT = -expm1(-GT);
    
    D  = -expm1(-2*GT);
    if sigma > 0
        D = D + expm1(-GT).^2 *sigma^2;
    end
    vec = ((x1-x0)/dt - (x0-xm1)/dt*EGT - fun(p, x0)/GM*M1EGT).^2;
    r =   A1*log(pi*B*D) + 1/B/D .* sum(vec, 'all');
end

function r = gsumfun_p(p, GM, x1, x0, xm1, fun, dt, kBT_m, sigma)
    A1    =  numel(x0)/2;
    B     = kBT_m*2; 
    GM    = abs(GM);
    GT    = GM*dt;
    M1EGT = -expm1(-GT);
    if GT < 1e-4
        D  = 2/3*GT^3 - 0.5*GT^4;
    else
        D = 2*GT + 4*expm1(-GT) - expm1(-2*GT);
    end
    if sigma > 0
        D = D + expm1(-GT).^2 *sigma^2;
    end
    D = D / GM^2;
    vec = (x1-x0 - (x1-xm1)/2/GT*M1EGT - fun(p, x0)/GM^2*(GT - M1EGT)).^2;
    r =   A1*log(pi*B*D) + 1/B/D .* sum(vec, 'all');
end


function r = gsumfun_c(p, GM, x1, x0, xm1, fun, dt, kBT_m, sigma)
    D0 = abs(kBT_m);
    DIMS = size(x0,2);

    GM = abs(GM);
    GT = GM*dt;
    EGT  = exp(-GT);
    EGT2 = EGT^2;
    M1EGT= -expm1(-GT);
    
    if GT < 1e-4
        D1  = 4/3 - 4/3*GT + 4/5*GT^2;
        D2  = 1/3 - 1/3*GT + 11/60*GT^2;
    else
        D1 = (2*GT-2  + (2*GT+2)*EGT2)/GT^3;
        D2 = (1-EGT2-2*GT*EGT)/GT^3;
    end
    
    D0 = D0 / GM^2 * GT^3;
    vec = x1 - x0 * (1+EGT) + xm1 * EGT - fun(p, x0)/GM^2 * GT*M1EGT;

    N = numel(vec);

    if any(sigma > 0)
        C = zeros(3, N);
        fder = 0;
        if numel(sigma) == 1
            sigma = sigma*ones(DIMS,1);
        end
        for kk = 1:DIMS
            deltax = zeros(1, DIMS);
            deltax(kk) = sigma(kk)/10;            
            fder = fder + (fun(p, x0 + deltax) - fun(p, x0 - deltax))/2/deltax(kk);            
        end
        fder = fder / GM^2 * GT*M1EGT;
        
        Q = 0;
        logdetCC = 0;
        for kk = 1:DIMS
            C(1,:)       = D1 + sigma(kk)^2/D0 * (2 + 2*EGT + 2*EGT2 + 2*(1+EGT)*fder(:,kk) + fder(:,kk).^2);        
            C(2,1:end-1) = D2 - sigma(kk)^2/D0 * ((1+EGT)^2 + fder(1:end-1, kk) + fder(2:end, kk));
            C(3,1:end-2) =      sigma(kk)^2/D0 * EGT;
            [QQ, ldC] = solveBeepsisArg(C, vec);

            logdetCC = logdetCC + ldC;
            Q = Q  + QQ;
        end
    else
        C = zeros(2, N);
        C(1,:)       = D1;        
        C(2,1:end-1) = D2;
        [Q, logdetCC] = solveBeepsisArg(C, vec);
    end
               
    r =   0.5*( N*log(D0*2*pi) + logdetCC + 1/D0 .* Q);
end


