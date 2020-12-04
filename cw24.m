clear

fre_equations_680029911
load('udBranch.mat')

T = 2.*pi./omega;
N = floor(T/0.03);

magDisp = 1e-5;

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

pd1 = MySolve(A,uIni,JA);
pd1 = pd1(1:3);

%4b
figure()
scatter(udBranch(3,:),udBranch(1,:),15,udStab,'filled')

pd = pd1;

periodJ = 1;

[eVecs2,~] = eigs(JM(pd(1:end-1),pd(end)),1,-1);


FA = @(u) PeriodNSys(u(1:end-1),u(end),M,2^periodJ);
JFA = @(u) MyJacobian(@(y) FA([y;u(end)]),u(1:end-1),1e-6);

uIni2 = [pd(1:end-1)+magDisp*eVecs2;pd(1:end-1)-magDisp*eVecs2;pd(end)];

p2List = MyTrackCurve(FA,uIni2,[eVecs2;-eVecs2;0],'stop',@(y) y(end) > 3,'sMax',1e-3);

i = 1;
closeEvecIndex = 1;
closeEvec = [];
closeEval = inf;
while i <= size(p2List,2)
    
    [pEVec,pEVal] = eigs(JFA(p2List(:,i)),1,-1);
    
    eValDist = abs(-1-pEVal);
    if(eValDist < abs(-1-closeEval))
        closeEvecIndex = i;
        closeEvec = pEVec;
        closeEval = pEVal;
    end
    
    i = i + 1;
end

PDBSys = @(u) [FA(u(1:length(uIni2)));...
    JFA(u(1:length(uIni2)))*u(length(uIni2)+1:end)+u(length(uIni2)+1:end);...
    u(length(uIni2)+1:end)'*u(length(uIni2)+1:end)-1];

pdbStartU = [p2List(:,closeEvecIndex);closeEvec];

pdb3 = MySolve(PDBSys,pdbStartU,@(y) MyJacobian(PDBSys,y,1e-6));

hold on
scatter(p2List(end,:),p2List(1,:),15)
scatter(p2List(end,:),p2List(3,:),15)
plot(pdb3(5),pdb3(1),'kx')
plot(pdb3(5),pdb3(3),'kx')
hold off

PeriodNGuSys(p2List(1:4,10),p2List(end,10),M,2,1)

%%Functions

function J = CompositeJacobian(u,a,f,n)

nU = length(u)/n;

J = eye(nU);

i = 1;
while i <= n
    
    J = J * MyJacobian(@(y) f(y,a),u(nU*(i-1)+1:i*nU),1e-6);  
    i = i + 1;
    
end

end

function fu = PeriodNSys(u,a,f,n)

nU = length(u)/n;

fu = nan(nU*n,1);

i = 1;
while i <= n
    fu(nU*(i-1)+1:i*nU) = f(u(nU*(i-1)+1:i*nU),a) - u(mod(i,n)*nU+1 : (mod(i,n)+1)*nU);
    i = i + 1;
end

end

function gu = PeriodNGuSys(u,a,f,n,s)

nU = length(u)/n;
gu = zeros(n*nU);

if(n == 1)
    gu = MyJacobian(f(u,a),u,1e-6);
else
    
    i = 1;
    while i <= n
        
        indexRange = nU*(i-1)+1:i*nU;
        gu(indexRange,indexRange) = MyJacobian(@(y) f(y,a),u(indexRange),1e-6);
        
        i = i + 1;
        
    end
    
    I = -eye((n-1)*nU);
    I = padarray(I',nU,0,'pre')';
    I = padarray(I,2,0,'post');
    
    gu = gu + I;
    
    gu(end-nU+1:end,1:nU) = s*eye(nU);
    
end

end

function ga = PeriodNGaSys(u,a,f,n)

nU = length(u)/n;

ga = nan(nU*n,1);

i = 1;
while i <= n
    
    ga(nU*(i-1)+1:i*nU) = MyJacobian(@(y) f(u(nU*(i-1)+1:i*nU),y),a,1e-6);
    
    i = i + 1;
    
end

end