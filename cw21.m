clear

fre_equations_680029911

%Equations defining fixed points at a = 0
fr = @(r) pi .^2 .* r.^4 - Jbar .* r.^3 - eta .* r.^2 - (((Delta + Gamma.*r).^2)./(4.*pi.^2));
fv = @(r) -(Delta + Gamma.*r)./(2.*pi.*r);

%Stroboscopic map
T = 2.*pi./omega;
N = floor(T/0.03);
M = @(u0,a) MyIVP(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');

%Initial guesses for r at fixed points a = 0
rIniGuesses = [0.086,0.449,1.042];
rIni = [0,0,0];

%Converge guesses
i = 1;
for ri = rIniGuesses
    rIni(i) = MySolve(fr,ri,@(r) MyJacobian(fr,r,1e-6));
    i = i + 1;
end

%Calculate v from r values
vIni = fv(rIni);

%Stable points at a = 0;
uIni = [rIni;vIni]; 

%Function for period 1 orbit
f = @(x) M(x(1:2),x(3)) - x(1:2);

%Track from points at a=0
yList = [];
i = 1;
while i <= size(uIni,2)
    yList(:,:,i) = MyTrackCurve(f,[uIni(:,i);0],[0;0;1],'stop',@(y) y(3) > 3,'sMax',5e-2);
    i = i + 1;
end


%%Compute stability; 0 = saddle, 1 = stable, 2 = unstable
stab = NaN(size(uIni,2),length(yList));

dM = @(x) MyJacobian(@(y) M(y,x(3)),x(1:2),1e-6);

for j = 1:size(uIni,2)
    
    for i = 1:length(yList(1,:,j))
        
        if(~isnan(yList(:,i,j)))
            
            %Calculate eigenvalues
            eVals = eigs(dM(yList(:,i,j)));
            
            %determine stability
            %1 Stable
            %2 Unstable
            %3 Saddle
            
            if(all(abs(eVals) < 1))
                stab(j,i) = 1;
            elseif(all(abs(eVals) > 1))
                stab(j,i) = 2;
            else
                stab(j,i) = 3;
            end
            
        end
        
    end
end

figure()

%Colours used in plot
cMap = [[46,140,41];[245,35,8];[204,204,25]]./255;
colormap(cMap)

scatter(yList(3,:,1),yList(1,:,1),15,stab(1,:),'filled')
hold on
scatter(yList(3,:,2),yList(1,:,2),15,stab(2,:),'filled')
scatter(yList(3,:,3),yList(1,:,3),15,stab(3,:),'filled')

%Create plot legend
leg = zeros(3, 1);
leg(1) = plot(NaN,NaN,'o','MarkerFaceColor',[46,140,41]./255,'MarkerEdgeColor',[46,140,41]./255);
leg(2) = plot(NaN,NaN,'o','MarkerFaceColor',[245,35,8]./255,'MarkerEdgeColor',[245,35,8]./255);
leg(3) = plot(NaN,NaN,'o','MarkerFaceColor',[204,204,25]./255,'MarkerEdgeColor',[204,204,25]./255);
legend(leg, 'Stable','Unstable','Saddle','Location','northwest');

hold off

xlim([0,3])
xlabel('a')
ylabel('r')
title('Stability of fixed points in a-r plane')

figure()

colormap(cMap)

plot(yList(3,:,1),yList(2,:,1));
scatter(yList(3,:,1),yList(2,:,1),15,stab(1,:),'filled')
hold on
scatter(yList(3,:,2),yList(2,:,2),15,stab(2,:),'filled')
scatter(yList(3,:,3),yList(2,:,3),15,stab(3,:),'filled')

%Create plot legend
leg = zeros(3, 1);
leg(1) = plot(NaN,NaN,'o','MarkerFaceColor',[46,140,41]./255,'MarkerEdgeColor',[46,140,41]./255);
leg(2) = plot(NaN,NaN,'o','MarkerFaceColor',[245,35,8]./255,'MarkerEdgeColor',[245,35,8]./255);
leg(3) = plot(NaN,NaN,'o','MarkerFaceColor',[204,204,25]./255,'MarkerEdgeColor',[204,204,25]./255);
legend(leg, 'Stable','Unstable','Saddle','Location','west');

hold off

xlim([0,3])

xlabel('a')
ylabel('v')
title('Stability of fixed points in a-v plane')

%%1c

%Number of iterations
k = 5000;

%Initial conditions
a1 = 2.5;
u1 = uIni(:,3) + 1e-2;

%Compute lyapunov exponents
[lambda,rDiag,x] = LyapunovQR(@(u) M(u,a1),u1,k);

%Compute lyapunov exponents after each iteration
lyapExp = cumsum(log(rDiag),2)./[1:k;1:k];

kRange = 100:k; %Range of iterates to use, neglecting transients

figure()
plot(kRange,x(1,kRange),'o',kRange,x(2,kRange),'x')
xlabel('Itearation')
ylabel('Value')
legend('r','v')
title(sprintf('r and v at a = %.4g',a1))

figure()
plot(x(1,kRange),x(2,kRange),'o')
xlabel('r')
ylabel('v')
title(sprintf('Trajectory for a = %.4g',a1))

figure()
plot(1:k,lyapExp(1,:),1:k,lyapExp(2,:))
xlabel('Iteration')
ylabel('Lambda')
title(sprintf('Lyapunov exponents over iterations for a = %.4g',a1))

disp("Lyupanov Exponents at 5000 iterations:")
disp(lambda)

udBranch = yList(:,:,3);
udStab = stab(3,:);
save('udBranch','udBranch','udStab')
