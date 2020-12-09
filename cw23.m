clear

fre_equations_680029911

load('basins.mat')

%Stroboscopic map, and inverse stroboscopic map
T = 2.*pi./omega;
N = floor(T/0.03);

F = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,T],N,'dp45');
Finv = @(u0,a) MyIVPVec(@(t,u) rhs(u,a,t),u0,[0,-T],N,'dp45');

mapLineIterates = 4;

aLoop = 1;
while aLoop <= length(aList)
    
    a = aList(aLoop);
    
    %Map and inverse, updated definition for mapline
    M = @(u0) F(u0,a);
    JM = @(u0) MyJacobian(M,u0,1e-6);
    Minv = @(u0) Finv(u0,a);
    JMinv = @(u0) MyJacobian(Minv,u0,1e-6);
    
    roots = [];
    
    %Load saved roots
    i = 1;
    while i <= size(savedRoots,2)
        
        if(savedRoots(3,i) == a)
            roots = [roots,savedRoots([1,2],i)];
        end
        
        i = i + 1;
        
    end
    
    fprintf("\na = %4g\n",a)
    
    %determine and print root stabilities, as in 1a
    i = 1;
    while i <= size(roots,2)
        
        fprintf("\nRoot u_%c\n",char(i+97))
        fprintf("r = %.4g, v = %.4g\n\n",roots(1,i),roots(2,i))
        
        [eVecs,eVals] = eigs(JM(roots(:,i)));
        
        disp("Eigenvectors:")
        disp(eVecs)
        disp("Eigenvalues:")
        disp(diag(eVals))
        
        %Properties of u_c
        if(i==2)
            eVecsUC = eVecs;
            eValsUC = diag(eVals);
        end
        %Properties of u_d
        if(i==3)
            unstableEVecD = eVecs(:,1);
        end
        
        %Determine stability
        if(all(abs(eVals) < 1))
            fprintf("Stable\n")
        elseif(all(abs(eVals) > 1))
            fprintf("Unstable\n")
        else
            fprintf("Saddle\n")
        end
        
        i = i + 1;
        
    end
    
    %Define stable and unstable directions
    unstableEVecC = eVecsUC(:,1);
    stableEVecC = eVecsUC(:,2);
    
    %Initialize line segment for stable and unstable manifolds
    sOldU = [0,1];
    xOldU = [roots(:,2)+1e-2*unstableEVecC,roots(:,2)-5e-3*unstableEVecC];
    sOldS = [0,1];
    xOldS = [roots(:,2)+1e-2*stableEVecC,roots(:,2)-5e-3*stableEVecC];
    
    %Develop and refine manifolds
    i = 1;
    while i <= mapLineIterates
        
        [yNewU,xNewU,sNewU] = MapLine(M,xOldU,sOldU,0.01,deg2rad(5));
        [yNewS,xNewS,sNewS] = MapLine(Minv,xOldS,sOldS,0.01,deg2rad(5));
        
        xOldU = yNewU;
        sOldU = sNewU;
        xOldS = yNewS;
        sOldS = sNewS;
        
        i = i + 1;
        
    end
        
    figure()
    
    %load basin image
    convRegion = savedBasins(:,:,aLoop);
    imagesc(r,v,convRegion)
    
    colormap([[0,0,0];[0,1,0]]);
    set(gca,'YDir','normal')
    
    hold on
    %Plot roots
    i = 1;
    while i <= size(roots,2)
        plot(roots(1,i),roots(2,i),'*w');
        i = i + 1;
    end
    %Plot manifolds
    plot(yNewU(1,:),yNewU(2,:),"r",'linewidth',1)
    plot(yNewS(1,:),yNewS(2,:),"b",'linewidth',1)
    %Compute and plot trajectory from u_d
    if(a ~= 0)
        i = 1;
        uList = nan(2,1000);
        uList(:,1) = roots(:,3)+1e-3*unstableEVecD;
        while i < 1000
            uList(:,i+1) = M(uList(:,i));
            i = i + 1;
        end
        plot(uList(1,:),uList(2,:),'.','MarkerSize',5)
    end
    hold off
    title(['a = ', num2str(a)])
    
    aLoop = aLoop + 1;
    
end