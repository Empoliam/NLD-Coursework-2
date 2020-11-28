function o = MyJacobian(f,x,h)
%MyJacobian Estimate the value of the jacobian df/dx using a centered finite difference
%   Input:
%   
%   f - Column vector, function handles
%   x - Column vector, numeric, point at which to evaluate
%   h - Numeric, Finite difference stepsize
%
%   Output:
%   
%   o - Jacobian estimate
%

kMax = length(x); %number of columns
nf = length(f(x));  %Number of rows
I = eye(kMax); %Identity matrix, used to select appropriate x(n) to differentiate by

o = zeros(nf,kMax); %Initialise output

%Column by column calculation
for k = 1:kMax
    o(:,k) = (f(x + I(:,k).*h) - f(x - I(:,k).*h))/(2.*h); %centered finite difference
end

end