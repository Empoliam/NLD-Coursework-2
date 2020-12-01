clear

fre_equations_680029911

T = 2.*pi./omega;
N = 50;

M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');

nIterates = 25;

rDivs = 400;
vDivs = 400;

rMin = -0.25;
rMax = 2.5;

vMin = -2.5;
vMax = 2.5;

convergedTolerance = 1e-2;

aList = [0,2.1,2.6];

for a = aList
        
    %Find roots
    
    f = @(u) M(u,a) - u;
    df = @(u) MyJacobian(f,u,1e-6);
    
    rCheck = linspace(rMin,rMax,10);
    vCheck = linspace(vMin,vMax,10);
    
    [RCheck,VCheck] = meshgrid(rCheck,vCheck);
    
    RCheck = reshape(RCheck,1,[]);
    VCheck = reshape(VCheck,1,[]);
    
    UCheck = [RCheck;VCheck];
    
    rootConverged = nan(2,size(UCheck,2));
    
    i = 1;
    while i <= size(UCheck,2)
        rootConverged(:,i) = MySolve(f,UCheck(:,i),df,'maxIter',20);
        i = i + 1;
    end
    
    rootConverged = rootConverged(:,all(~isnan(rootConverged)));
    roots = uniquetol(rootConverged',1e-4,'ByRows',true)';
    
    %find basin
    
    r = linspace(rMin,rMax,rDivs);
    v = linspace(vMin,vMax,vDivs);
    
    [R,V] = meshgrid(r,v);
    
    R = reshape(R,1,[]);
    V = reshape(V,1,[]);
    
    U = [R;V];
    
    tic
    i = 1;
    while i <= nIterates
        [uend, t, ut] = M(U,a);
        U = uend;
        i = i + 1;
    end
    toc
    
    convRegion = zeros(1,size(uend,2));
    convDist = zeros(1,size(uend,2));
    i = 1;
    while i <= size(uend,2)
        
        if(norm(uend(:,i)-roots(:,1)) < convergedTolerance)
            convRegion(i) = 1;
        end
        
        convDist(i) = norm(uend(:,i)-roots(:,1));
        
        i = i + 1;
    end
    
    convRegion = reshape(convRegion,length(v),length(r));
    convDist = reshape(convDist,length(v),length(r));
    
    trajectoryUC = nan(2,N);
    trajectoryUC(:,1) = roots(:,2) - 10e-3;
    
    trajectoryUD = nan(2,N);
    trajectoryUD(:,1) = roots(:,3) + 10e-3;
    
    i = 1;
    while(i < nIterates)
        trajectoryUC(:,i+1) = M(trajectoryUC(:,i),a);
        trajectoryUD(:,i+1) = M(trajectoryUD(:,i),a);
        i = i + 1;
    end
    
    figure()
    colormap('hot')
    imagesc(r,v,convDist)
    set(gca,'YDir','normal')
    hold on
    i = 1;
    while i <= size(roots,2)
        plot(roots(1,i),roots(2,i),'xw');
        i = i + 1;
    end
    hold off
    title(['a = ', num2str(a)])
    
    figure()
    imagesc(r,v,convRegion)
    colormap([[0,0,0];[0,1,0]]);
    set(gca,'YDir','normal')
    hold on
    i = 1;
    while i <= size(roots,2)
        plot(roots(1,i),roots(2,i),'xw');
        i = i + 1;
    end
    plot(trajectoryUC(1,:),trajectoryUC(2,:),'bo')
    plot(trajectoryUD(1,:),trajectoryUD(2,:),'ro')
    hold off
    title(['a = ', num2str(a)])
         
end