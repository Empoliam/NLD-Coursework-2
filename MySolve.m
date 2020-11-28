function [x,converged,J]=MySolve(f,x0,df,varargin)
%MySolve Estimate the root of a given system of equations using the Newton method
%   Input:
%   
%   f - Column vector, function handles, function to solve
%   x0 - Column vector, numeric, initial guess
%   df - Function handle, function returning jacobian of f
%   varargin - Optional inputs
%   Output:
%   
%   x - Root estimate
%   converged - Numeric, 0 for failed convergence, 1 for sucessful convergence 
%   J - Numeric matrix, Last computed jacobian
%

%Varargs parsing
parser = inputParser;

addRequired(parser,'f',@(func) isa(func,'function_handle'));
addRequired(parser,'x0',@isnumeric);
addRequired(parser,'df',@(func) isa(func,'function_handle'));
addOptional(parser,'maxIter',10,@isnumeric);                        %Maximum solver iterations
addOptional(parser,'tol',1e-8,@isnumeric);                          %Tolerance for root
parse(parser,f,x0,df,varargin{:})

maxIter = parser.Results.maxIter;
tol = parser.Results.tol;

%Initialise iteration
xn = x0;

iterations = 0;
converged = 0;

while (converged == 0) && (iterations < maxIter)
        
    x = xn;
    J = df(x);  %Compute jacobian at current guess
    
    xn = x - (J^(-1)) * f(x);   %Newton method
    
    if (all(abs(x-xn) < tol) || (norm(f(x)) < tol)) %Break loop once root is within tolerance
        converged = 1;
        x = xn;
        return
    end
    
    iterations = iterations + 1;
    
end

end