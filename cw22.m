clear

fre_equations_680029911

T = 2.*pi./omega;
N = 100;

M = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');

rDivs = 500;
vDivs = 500;

rMin = -0.25;
rMax = 2.5;

vMin = -2.5;
vMax = 2.5;

convergedTolerance = 1.25;

aList = [0,2.1,2.5];

for a = aList
    
    tic
    
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
    
    [uend, t, ut] = M(U,a);  
    
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
        
    figure()
    colormap('turbo')
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
    hold off
    title(['a = ', num2str(a)])
    
    toc
    
end