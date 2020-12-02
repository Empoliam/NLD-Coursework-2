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