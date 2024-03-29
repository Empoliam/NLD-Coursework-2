clear

fre_equations_680029911

T = 2.*pi./omega;
N = floor(T/0.03);

M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');

%Number of iterations to perform at each sample point
nIterates = 50;

%Resolution of r and v axes
rDivs = 500;
vDivs = 500;

%Range of axes
rMin = -0.25;
rMax = 2.5;

vMin = -2.5;
vMax = 2.5;

%Tolerance for a converged point
convergedTolerance = 1e-2;

%List of a values to use
aList = [0,2.1,2.6];
savedRoots = [];
savedBasins = [];

%Sample grid for locating fixed points
rCheck = linspace(rMin,rMax,10);
vCheck = linspace(vMin,vMax,10);

[RCheck,VCheck] = meshgrid(rCheck,vCheck);

RCheck = reshape(RCheck,1,[]);
VCheck = reshape(VCheck,1,[]);
UCheck = [RCheck;VCheck];

%Sample grid for checking convergence
r = linspace(rMin,rMax,rDivs);
v = linspace(vMin,vMax,vDivs);

[R,V] = meshgrid(r,v);

R = reshape(R,1,[]);
V = reshape(V,1,[]);

aLoop = 1;
while aLoop <= length(aList)
    
    a = aList(aLoop);
    
    %Find roots
    
    %Defining system for fixed points
    f = @(u) M(u,a) - u;
    df = @(u) MyJacobian(f,u,1e-6);
    
    rootConverged = nan(2,size(UCheck,2));
    
    %Solve for fixed points at each sample point
    i = 1;
    while i <= size(UCheck,2)
        rootConverged(:,i) = MySolve(f,UCheck(:,i),df,'maxIter',20);
        i = i + 1;
    end
    
    %Remove nans
    rootConverged = rootConverged(:,all(~isnan(rootConverged)));
    %Identify roots
    roots = uniquetol(rootConverged',1e-4,'ByRows',true)';
    
    %Save roots for future use
    savedRoots = [savedRoots,[roots;repelem(a,size(roots,2))]];
    
    %find basin
    
    %Create sample grid
    U = [R;V];

    %Iterate points on sample grid and save final location
    i = 1;
    while i <= nIterates
        [uend, t, ut] = M(U,a);
        U = uend;
        i = i + 1;
    end
        
    convRegion = zeros(1,size(uend,2));
    i = 1;
    while i <= size(uend,2)
        
        %Check if point is within tolerance of u_b
        if(norm(uend(:,i)-roots(:,1)) < convergedTolerance)
            convRegion(i) = 1;
        end
        
        i = i + 1;
    end
    
    convRegion = reshape(convRegion,length(v),length(r));
    
    %Save basin for future use
    savedBasins(:,:,aLoop) = convRegion;
    
    %Initialise trajectories from fixed points
    trajectoryUC = nan(2,N);
    trajectoryUC(:,1) = roots(:,2) - 10e-3;
    
    trajectoryUD = nan(2,N);
    trajectoryUD(:,1) = roots(:,3) + 10e-3;
    
    %Compute trajectories
    i = 1;
    while(i < 50)
        trajectoryUC(:,i+1) = M(trajectoryUC(:,i),a);
        trajectoryUD(:,i+1) = M(trajectoryUD(:,i),a);
        i = i + 1;
    end
    
    %Draw results
    figure()
    imagesc(r,v,convRegion)
    colormap([[0,0,0];[0,0.9,0]]);
    set(gca,'YDir','normal')
    hold on
    plot(trajectoryUC(1,:),trajectoryUC(2,:),'bo')
    plot(trajectoryUD(1,:),trajectoryUD(2,:),'ro')
    i = 1;
    while i <= size(roots,2)
        plot(roots(1,i),roots(2,i),'xw');
        i = i + 1;
    end
    hold off
    title(sprintf('Basin of attraction of ub for a = %.4g',a ))
    
    %print root locations
    i = 1;
    while i <= size(roots,2)
        fprintf("\nRoot u_%c\n",char(i+97))
        fprintf("r = %.4g, v = %.4g\n\n",roots(1,i),roots(2,i))
        
        i = i + 1;
    end
    
    aLoop = aLoop + 1;
    
end

save("basins",'aList','savedRoots','savedBasins','r','v')
