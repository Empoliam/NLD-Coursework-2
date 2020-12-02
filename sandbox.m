clear

x = -1:0.25:1;
y = x.^2;

yList = [x;y]

i = 5;

segA = yList(:,i) - yList(:,i-1);
segB = yList(:,i+1) - yList(:,i-1);
segAngle = acos(dot(segA./norm(segA),segB./norm(segB)));

rad2deg(segAngle)

plot(x,y)
