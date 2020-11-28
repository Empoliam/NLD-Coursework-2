function [lambda,rDiag,x]=LyapunovQR(M,xini,N)
%LYAPUNOVQR Compute the Lyapunov Exponents of the given map along a trajectory starting at xini

x = nan(length(xini),N);
x(:,1) = xini;

Q = eye(length(xini),N);

rDiag = nan(length(xini),N);

i = 1;
while i <= N
    
    A = MyJacobian(M,x(:,i),1e-8);
    [Q,R] = qr(A*Q);
    
    rSign = diag(sign(diag(R)));
    Q = Q*rSign;
    R = rSign*R;
    
    rDiag(:,i) = diag(R);
    
    if (i < N) 
        x(:,i+1) = M(x(:,i)) ;
    end
    
    i = i + 1;
    
end

lambda = sum(log(rDiag),2) ./ N;

end

