function [lambda,rDiag,x]=LyapunovQR(M,xini,N)
%LYAPUNOVQR Compute the Lyapunov Exponents of the given map along a trajectory starting at xini
%
% Input:
%
% M - Map function
% xini - initial condition for iteration
% N - Number of iterates
%
% Output:
% x - Iterated trajectory
% rDiag - Diagonals of R at each iteration
% lambda - Lyapunov exponents at iteration N

x = nan(length(xini),N);
x(:,1) = xini;

Q = eye(length(xini),N);

rDiag = nan(length(xini),N);

i = 1;
while i <= N
    
    %Calculate jacobian at x(i)
    A = MyJacobian(M,x(:,i),1e-8);
    %Decompose
    [Q,R] = qr(A*Q);
    
    %Force positive diagonal of QR decomposition
    rSign = diag(sign(diag(R)));
    Q = Q*rSign;
    R = rSign*R;
    
    %Store diagonal
    rDiag(:,i) = diag(R);
    
    %Calculate next iteration
    if (i < N) 
        x(:,i+1) = M(x(:,i)) ;
    end
    
    i = i + 1;
    
end

%Calculate final lyapunov exponent
lambda = sum(log(rDiag),2) ./ N;

end

