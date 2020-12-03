clear

fre_equations_680029911
load('udBranch.mat')

T = 2.*pi./omega;
N = floor(T/0.03);

F = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');
JF = @(u0,a) MyJacobian(@(u) F(u,a),u0,1e-6);

%4a
A = @(u) [F(u(1:2),u(3))-u(1:2);...
    JF(u(1:2),u(3))*u(4:5)+u(4:5);...
    u(4:5)'*u(4:5)-1];
JA = @(u) MyJacobian(A,u,1e-6);

%Locate stable point close to pd
closeIndex = find(udStab==3,1);
closeU = udBranch(:,closeIndex);
[closeEVec,~] = eigs(JF(closeU(1:2),closeU(3)),1,-1);

uIni = [closeU;closeEVec];

MySolve(A,uIni,JA)

%4b

[X,Y] = PeriodNSys([0.087099592220115;-1.941260383954402],0.01,F,5);

function [u,fu] = PeriodNSys(ui,a,f,n)

    nU = length(ui);
    
    u = nan(nU*n + 1,1);
    u(1:nU) = ui;
    
    fu = nan(nU*n,1);
    
    i = 1;
    while i < n
        
        u(nU*(i)+1:(i+1)*nU) = f(u(nU*(i-1)+1:(i)*nU),a);
                
        i = i + 1;
    end

    fu(nU*(n-1)+1:nU*n) = f(u(nU*(n-1)+1:nU*n),a) - u(1:nU);
    u(end) = a;
    
end