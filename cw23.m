clear

fre_equations_680029911

T = 2.*pi./omega;
N = 50;

M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');