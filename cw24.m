clear

fre_equations_680029911
load('udBranch.mat')

%Resolution of various functions
magJac = 1e-10;
magDisp = 1e-10;
magH = 1e-10;

%Maximum number of bifurcations to compute (up to 2^(periodMax+1))
periodMax = 7;

%Stroboscopic Map
T = 2.*pi./omega;
N = floor(T/0.03);
M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');
JM = @(u0,a) MyJacobian(@(u) M(u,a),u0,magJac);

%List of a values of bifurcations
apdList = zeros(1,periodMax+1);

%4a - First bifurcation

%System for bifurcation point
A = @(u) [M(u(1:2),u(3))-u(1:2);...
    JM(u(1:2),u(3))*u(4:5)+u(4:5);...
    u(4:5)'*u(4:5)-1];
JA = @(u) MyJacobian(A,u,magJac);

%Locate stable point close to pd
closeIndex = find(udStab==3,1); %Change in stability
closeU = udBranch(:,closeIndex);
[closeEVec,~] = eigs(JM(closeU(1:2),closeU(3)),1,-1);

%Solve for point
uIni = [closeU;closeEVec];
pd1 = MySolve(A,uIni,JA);
pd1 = pd1(1:3);

apdList(1) = pd1(end);

figure()

plot(udBranch(3,closeIndex-10:closeIndex),udBranch(1,closeIndex-10:closeIndex))
hold on
plot(pd1(3),pd1(1),'kx')
hold off

xlabel('a')
ylabel('r')
title('Bifurcation diagram')

drawnow

%4b
%Track p2 orbits
%Find evec close to -1
[eVecs2,~] = eigs(JM(pd1(1:end-1),pd1(end)),1,-1);

%System for p2 orbits
F2A = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^1);
JF2A = @(u) MyJacobian(@(y) F2A([y;u(end)]),u(1:end-1),magJac);

%Track p2 orbits
uIni2 = [pd1(1:end-1)+magDisp*eVecs2;pd1(1:end-1)-magDisp*eVecs2;pd1(end)];
p2List = MyTrackCurve(F2A,uIni2,[eVecs2;-eVecs2;0],'stop',@(y) y(end) > 3,'sMax',1e-2,'nMax',100);

branchList = p2List;

%Compute further orbits
for periodI = 1:periodMax
     
    tic
    
    %System for p-n orbit
    FA = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^periodI);
    
    %Jacobian estimates
    GU = @(u,s) PeriodNGuSys(u(1:end-1),u(end),M,2^periodI,s);
    GA = @(u) PeriodNGaSys(u(1:end-1),u(end),M,2^periodI);
    
    %Second order jacobian estimates
    GUU = @(y,z,h) PeriodNGuuSys(y(1:end-1),y(end),M,2^periodI,z,h);
    GUA = @(y,h) PeriodNGuaSys(y(1:end-1),y(end),M,2^periodI,h);
        
    %Locate point close to PD
    i = 1;
    closeIndex = 1;
    oldEvals = zeros(2,1);
    while i <= size(branchList,2)
        
        if(~isnan(branchList(:,i)))
            
            %Get eigenvalues at current point
            newEvals = eigs(CompositeJacobian(branchList(1:end-1,i),branchList(end,i),M,2^(periodI)));
            
            %Check if any eigs cross -1
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
    
    %Close u value
    closeU = branchList(:,closeIndex);
    
    %range of indices defining [y,a]
    yIndices = 1:length(closeU);
    
    %System for period doubling bifurcation
    PDBSys = @(u) [FA(u(yIndices));...
        GU(u(yIndices),1)*u(length(closeU)+1:end) ;...
        u(length(closeU)+1:end).'*u(length(closeU)+1:end)-1];
    
    %Initial guess for PD bifurcation location
    [zGuess,~] = eigs(GU(closeU,1),1,0);
    uGuess = [closeU;real(zGuess)];
    
    %Assemble matrix for jacobian estimate 
    JUL = @(u) GU(u(yIndices),-1);
    JUC = @(u) GA(u(yIndices));
    JUR = @(u) zeros(length(u(yIndices))-1);
    JCL = @(u) GUU(u(yIndices),u(length(closeU)+1:end),magH);
    JCC = @(u) GUA(u(yIndices),magH) * u(length(closeU)+1:end);
    JCR = @(u) GU(u(yIndices),1);
    JLL = @(u) zeros(1,length(closeU)-1);
    JLC = @(u) 0;
    JLR = @(u) 2*u(length(closeU)+1:end).';
    
    ApproxJ = @(u) [JUL(u),JUC(u),JUR(u);JCL(u),JCC(u),JCR(u);JLL(u),JLC(u),JLR(u)];
        
    %Solve system
    pdbSolve = MySolve(PDBSys,uGuess,ApproxJ,'maxIter',25,'tol',1e-8);
    pdby = pdbSolve(1:length(closeU));
    pdbz = pdbSolve(length(closeU)+1:end);
     
    %Add a value to list
    apdList(periodI+1) = pdby(end);
    
    %Plot old branch and new bifurcation
    hold on
    sizeN = (length(closeU)-1)/2^periodI;
    i = 1;
    while i < length(closeU)
        
        plot(branchList(end,1:closeIndex),branchList(i,1:closeIndex))
        plot(pdby(end),pdby(i),'kx')
        
        i = i + sizeN;
    end
    hold off
    drawnow
    
    %Track new orbits
 
    %System for periodI+1 period orbits
    FB = @(u) PeriodNFSys(u(1:end-1),u(end),M,2^(periodI+1));
    FBJEst = @(u) [PeriodNGuSys(u(1:end-1),u(end),M,2^(periodI+1),-1),PeriodNGaSys(u(1:end-1),u(end),M,2^(periodI+1))];
    
    %initial guess
    uIniB = [pdby(1:end-1) + magDisp *pdbz; pdby(1:end-1) - magDisp * pdbz; pdby(end)];
    tanIniB = [pdbz;-pdbz;0];
    
    %Track new orbits
    branchList = MyTrackCurve(FB,uIniB,tanIniB,'stop',@(y) y(end) > 3,'sMax',(1e-2)/(2^(periodI)),'nMax',40,'userdf',FBJEst,'stepsize',1e-4);
    
    toc
    
end

%Plot final orbits
periodI = periodI + 1;
hold on
i = 1;
while i < sizeN*2^periodI
    
    plot(branchList(end,:),branchList(i,:))
    
    i = i + sizeN;
end
hold off

%Calculate ratios
i = 1;
while i < periodMax
    aRatio = (apdList(i+1)-apdList(i))/(apdList(i+2)-apdList(i+1));
    disp(aRatio)
    i = i + 1;
end

%%Functions

%Jacobian Jn*.....*j2*j1
function J = CompositeJacobian(u,a,f,n)

nU = length(u)/n;

J = eye(nU);

i = 1;
while i <= n
    
    J = MyJacobian(@(y) f(y,a),u(nU*(i-1)+1:i*nU),1e-10) * J;
    i = i + 1;
    
end

end

%System F in notes
function fu = PeriodNFSys(u,a,f,n)

nU = length(u)/n;

fu = nan(nU*n,1);

i = 1;
while i <= n
    fu(nU*(i-1)+1:i*nU) = f(u(nU*(i-1)+1:i*nU),a) - u(mod(i,n)*nU+1 : (mod(i,n)+1)*nU);
    i = i + 1;
end

end

%System Gu in notes
function gu = PeriodNGuSys(u,a,f,n,s)

nU = length(u)/n;
gu = zeros(n*nU);

i = 1;
while i <= n
    
    %Jacobians on diagonal
    indexRange = nU*(i-1)+1:i*nU;
    gu(indexRange,indexRange) = MyJacobian(@(y) f(y,a),u(indexRange),1e-10);
    
    i = i + 1;
    
end

%Generate off diagonal identity matrices
I = -eye((n-1)*nU);
I = padarray(I',nU,0,'pre')';
I = padarray(I,nU,0,'post');

gu = gu + I;

%Bottom left corner
gu(end-nU+1:end,1:nU) = s*eye(nU);

end

%System Ga in notes
function ga = PeriodNGaSys(u,a,f,n)

nU = length(u)/n;

ga = nan(nU*n,1);

i = 1;
while i <= n
    
    ga(nU*(i-1)+1:i*nU) = MyJacobian(@(y) f(u(nU*(i-1)+1:i*nU),y),a,1e-10);
    
    i = i + 1;
    
end

end

%System Guu in notes, centre difference estimate for d/dy Gu z
function guu = PeriodNGuuSys(u,a,f,n,z,h)

A = PeriodNGuSys(u+h*z,a,f,n,1);
B = PeriodNGuSys(u-h*z,a,f,n,1);

guu = (1/(2*h)) * (A-B);

end

%System Gua in notes, centre difference estimate for d/dy Gu z
function gua = PeriodNGuaSys(u,a,f,n,h)

A = PeriodNGuSys(u,a+h,f,n,1);
B = PeriodNGuSys(u,a-h,f,n,1);

gua = (1/(2*h)) * (A-B);

end