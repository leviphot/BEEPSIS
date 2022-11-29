function hessian = numhessian(fun,x)
% function Hessian = numhessian(fun,x)
%
% numerically calculates Hessian of the function fun in point x
% addaptd from sub-function finDiffHessian in fminsub
%
% INPUTS
%    fun    function handle
%           fun takes one vector input and returns one number
%    x      vector (point in N dim space) where hessian iss calulated
% 
% OUTPUT
%  hessian  nuemrically evaluated Hessian matrix


numberOfVariables = numel(x);
hessian = zeros(numberOfVariables);
f = fun(x);

% Define stepsize
CHG = eps^(1/4)*sign(x).*max(abs(x),1);


% Calculate the upper triangle of the finite difference Hessian element
% by element, using only function values. The forward difference formula
% we use is
%
% Hessian(i,j) = 1/(h(i)*h(j)) * [f(x+h(i)*ei+h(j)*ej) - f(x+h(i)*ei)
%                          - f(x+h(j)*ej) + f(x)]                   (2)
%
% The 3rd term in (2) is common within each column of Hessian and thus
% can be reused. We first calculate that term for each column and store
% it in the row vector fplus_array.
fplus_array = zeros(1,numberOfVariables);
for j = 1:numberOfVariables
    xplus = x;
    xplus(j) = x(j) + CHG(j);
    % evaluate
    fplus = fun(xplus);
    fplus_array(j) = fplus;
end

for i = 1:numberOfVariables
    % For each row, calculate the 2nd term in (2). This term is common to
    % the whole row and thus it can be reused within the current row
    xplus = x;
    xplus(i) = x(i) + CHG(i);
    % evaluate
    fplus_i = fun(xplus);
    
    for j = i:numberOfVariables   % start from i: only upper triangle
        % Calculate the 1st term in (2); this term is unique for each element
        % of Hessian and thus it cannot be reused.
        xplus = x;
        xplus(i) = x(i) + CHG(i);
        xplus(j) = xplus(j) + CHG(j);
        % evaluate
        fplus = fun(xplus);
        hessian(i,j) = (fplus - fplus_i - fplus_array(j) + f)/(CHG(i)*CHG(j));
    end
end

% Fill in the lower triangle of the Hessian
hessian = hessian + triu(hessian,1)';
