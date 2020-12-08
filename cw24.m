clear

fre_equations_680029911
load('udBranch.mat')

T = 2.*pi./omega;
N = floor(T/0.03);

magDisp = 1e-8;
magH = 1e-8;

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
%Track p2 orbits
figure()

[eVecs2,~] = eigs(JM(pd1(1:end-1),pd1(end)),1,-1);
F2A = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^1);
JF2A = @(u) MyJacobian(@(y) F2A([y;u(end)]),u(1:end-1),1e-8);
uIni2 = [pd1(1:end-1)+magDisp*eVecs2;pd1(1:end-1)-magDisp*eVecs2;pd1(end)];
p2List = MyTrackCurve(F2A,uIni2,[eVecs2;-eVecs2;0],'stop',@(y) y(end) > 3,'sMax',1e-2,'nMax',100);
% hold on
% scatter(p2List(end,:),p2List(1,:),15)
% scatter(p2List(end,:),p2List(3,:),15)
% hold off

branchList = p2List;

%Compute further orbits
for periodI = 1:4
    tic
   
    periodI
    
    FA = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^periodI);
    GU = @(u,s) PeriodNGuSys(u(1:end-1),u(end),M,2^periodI,s);
    GA = @(u) PeriodNGaSys(u(1:end-1),u(end),M,2^periodI);
    GUU = @(y,z,h) PeriodNGuuSys(y(1:end-1),y(end),M,2^periodI,z,h);
    GUA = @(y,h) PeriodNGuaSys(y(1:end-1),y(end),M,2^periodI,h);
    
    i = 1;
    closeIndex = 1;
    oldEvals = zeros(2,1);
    while i <= size(branchList,2)
        
        if(~isnan(branchList(:,i)))
            
            newEvals = eigs(CompositeJacobian(branchList(1:end-1,i),branchList(end,i),M,2^(periodI)));
            
            crossingCriteria = any(newEvals(oldEvals <= -1) >= -1) || any(newEvals(oldEvals >= -1) <= -1);
            if(crossingCriteria)
                               
                closeIndex = i;
                break
            else
                oldEvals = newEvals;
            end
            
            
        end
        i = i + 1;
    end
    
    closeU = branchList(:,closeIndex);
    
    yIndices = 1:length(closeU);
    
    PDBSys = @(u) [FA(u(yIndices));...
        GU(u(yIndices),1)*u(length(closeU)+1:end) ;...
        u(length(closeU)+1:end)'*u(length(closeU)+1:end)-1];
    
    [zGuess,~] = eigs(PeriodNGuSys(closeU(1:end-1),closeU(end),M,2^(periodI),1),1,'smallestreal');
    uGuess = [closeU;zGuess];
    
    JUL = @(u) GU(u(yIndices),-1);
    JUC = @(u) GA(u(yIndices));
    JUR = @(u) zeros(length(u(yIndices))-1);
    JCL = @(u) GUU(u(yIndices),u(length(closeU)+1:end),magH);
    JCC = @(u) GUA(u(yIndices),magH) * u(length(closeU)+1:end);
    JCR = @(u) GU(u(yIndices),1);
    JLL = @(u) zeros(1,length(closeU)-1);
    JLC = @(u) 0;
    JLR = @(u) 2*u(length(closeU)+1:end)';
    
    ApproxJ = @(u) [JUL(u),JUC(u),JUR(u);JCL(u),JCC(u),JCR(u);JLL(u),JLC(u),JLR(u)];
     
%     imag(ApproxJ(uGuess))
%     imag(MyJacobian(PDBSys,uGuess,1e-10))
    
    pdbSolve = MySolve(PDBSys,uGuess,@(u) MyJacobian(PDBSys,u,1e-10));
    pdby = pdbSolve(1:length(closeU));
    pdbz = pdbSolve(length(closeU)+1:end);
    
    hold on
    scatter(branchList(end,:),branchList(1,:),15)
%     scatter(branchList(end,:),branchList(3,:),15)
    plot(pdby(end),pdby(1),'kx')
%     plot(pdb3(length(closeU)),pdb3(3),'kx')
    hold off
    drawnow
    
    %Track new orbits
    
    FB = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^(periodI+1));
    
    uIniB = [pdby(1:end-1) + magDisp *pdbz; pdby(1:end-1) - magDisp * pdbz; pdby(end)];
    tanIniB = [pdbz;-pdbz;0];
    
    branchList = MyTrackCurve(FB,uIniB,tanIniB,'stop',@(y) y(end) > 3,'sMax',1e-2,'nMax',25);
    
    toc
    
end

%%Functions

function J = CompositeJacobian(u,a,f,n)

nU = length(u)/n;

J = eye(nU);

i = 1;
while i <= n
    
    J = J * MyJacobian(@(y) f(y,a),u(nU*(i-1)+1:i*nU),1e-8);
    i = i + 1;
    
end

end

function fu = PeriodNFSys(u,a,f,n)

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



i = 1;
while i <= n
    
    indexRange = nU*(i-1)+1:i*nU;
    gu(indexRange,indexRange) = MyJacobian(@(y) f(y,a),u(indexRange),1e-8);
    
    i = i + 1;
    
end

I = -eye((n-1)*nU);
I = padarray(I',nU,0,'pre')';
I = padarray(I,nU,0,'post');

gu = gu + I;

gu(end-nU+1:end,1:nU) = s*eye(nU);


end

function ga = PeriodNGaSys(u,a,f,n)

nU = length(u)/n;

ga = nan(nU*n,1);

i = 1;
while i <= n
    
    ga(nU*(i-1)+1:i*nU) = MyJacobian(@(y) f(u(nU*(i-1)+1:i*nU),y),a,1e-10);
    
    i = i + 1;
    
end

end

function guu = PeriodNGuuSys(u,a,f,n,z,h)

A = PeriodNGuSys(u+h*z,a,f,n,1);
B = PeriodNGuSys(u-h*z,a,f,n,1);

guu = (1/(2*h)) * (A-B);

end

function gua = PeriodNGuaSys(u,a,f,n,h)

A = PeriodNGuSys(u,a+h,f,n,1);
B = PeriodNGuSys(u,a-h,f,n,1);

gua = (1/(2*h)) * (A-B);

end