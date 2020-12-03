clear

fre_equations_680029911
load('udBranch.mat')

T = 2.*pi./omega;
N = floor(T/0.03);

M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');
JM = @(u0,a) MyJacobian(@(u) M(u,a),u0,1e-6);

%4a
A = @(u) [M(u(1:2),u(3))-u(1:2);...
    JM(u(1:2),u(3))*u(4:5)+u(4:5);...
    u(4:5)'*u(4:5)-1];
JA = @(u) MyJacobian(A,u,1e-6);

%Locate stable point close to pd
closeIndex = find(udStab==3,1);
closeU = udBranch(:,closeIndex);
[closeEVec,~] = eigs(JM(closeU(1:2),closeU(3)),1,-1);

uIni = [closeU;closeEVec];

apd1 = MySolve(A,uIni,JA);
apd1 = apd1(1:3);

%4b
magDisp = 1e-5;

figure()
scatter(udBranch(3,:),udBranch(1,:),15,udStab,'filled')

[eVecs2,~] = eigs(JM(apd1(1:2),apd1(3)),1,-1);

periodJ = 1;
FA = @(u) PeriodNSys(u(1:end-1),u(end),M,2^periodJ);

uIni2 = [apd1(1:end-1)+magDisp*eVecs2;apd1(1:end-1)-magDisp*eVecs2;apd1(end)];

p2List = MyTrackCurve(FA,uIni2,[eVecs2;-eVecs2;0],'stop',@(y) y(end) > 3,'sMax',1e-2);

M(M(p2List(1:2,10),p2List(end,10)),p2List(end,10)) - p2List(1:2,10)

hold on
scatter(p2List(end,:),p2List(1,:),15)
scatter(p2List(end,:),p2List(3,:),15)
hold off

function fu = PeriodNSys(u,a,f,n)

    nU = length(u)/n;
        
    fu = nan(nU*n,1);
    
    i = 1;
    while i <= n
        fu(nU*(i-1)+1:i*nU) = f(u(nU*(i-1)+1:i*nU),a) - u(mod(i,n)*nU+1 : (mod(i,n)+1)*nU);
        i = i + 1;
    end
 
end