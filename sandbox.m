clear


f = @(x,t) 2.*x;

x0 = 1:10;

[xend, t, xt] = MyIVPVec(f,x0,[0,1],100,'dp45');