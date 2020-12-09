function [yNew,xNew,sNew] = MapLine(M,x,s,maxDist,maxAngle)
%MAPLINE Apply the given map to the input line segment, then refine to reduce complexity
%
%Input:
% M - Map function
% x - input segment
% s - Array s such that x(s(i)) = x(i)
% maxDist - Maximum distance between succesive points
% maxAngle - Maximum angle formed by two successive line segments
%
% Output:
% yNew - refined mapping of input x
% xNew - refined input x
% sNew - refined input s

%Initialize
yNew = M(x);
sNew = s;
xNew = x;

iMax = length(sNew);
i = 1;
while i < iMax
    
    %Insert new point if distance between points is too great
    while( norm(yNew(:,i+1) - yNew(:,i)) > maxDist)
        
        sStar = (sNew(i) + sNew(i+1)) ./ 2;
        xStar = interp1(sNew,xNew',sStar,'spline')';
        yStar = M(xStar);
        
        sNew = ColumnInsert(sStar,sNew,i+1);
        xNew = ColumnInsert(xStar,xNew,i+1);
        yNew = ColumnInsert(yStar,yNew,i+1);
        
    end
    
    if(i ~= 1)
         
        %calculate angle between points
        segA = yNew(:,i) - yNew(:,i-1);
        segB = yNew(:,i+1) - yNew(:,i-1);
        segAngle = acos(dot(segA./norm(segA),segB./norm(segB)));
        
        %Insert new points to left and right if criteria not met
        if(segAngle > maxAngle)
            
            sStarPlus = (sNew(i) + sNew(i + 1)) ./ 2;
            xStarPlus = interp1(sNew,xNew',sStarPlus,'spline')';
            yStarPlus = M(xStarPlus);
            sStarMinus = (sNew(i) + sNew(i - 1)) ./ 2;
            xStarMinus = interp1(sNew,xNew',sStarMinus,'spline')';
            yStarMinus = M(xStarMinus);
            
            ColumnInsert(sStarPlus,sNew,i+1);
            ColumnInsert(xStarPlus,xNew,i+1);
            ColumnInsert(yStarPlus,yNew,i+1);
            ColumnInsert(sStarMinus,sNew,i);
            ColumnInsert(xStarMinus,xNew,i);
            ColumnInsert(yStarMinus,yNew,i);
            
        end
        
    end
    
    iMax = length(sNew);
    i = i + 1;
    
end

iMax = length(sNew);
i = 2;
while i < iMax
    
    %Calculate angle between points
    segA = yNew(:,i) - yNew(:,i-1);
    segB = yNew(:,i+1) - yNew(:,i-1);
    segAngle = acos(dot(segA./norm(segA),segB./norm(segB)));
   
    %Check if angle is within criteria
    if(segAngle < maxAngle)
        %check if distance between adjacent points fits distance criteria
        if(norm(yNew(:,i+1) - yNew(:,i-1)) < maxDist)
            %Remove point if simplified segment fits both criteria
            sNew = ColumnRemove(sNew,i);
            xNew = ColumnRemove(xNew,i);
            yNew = ColumnRemove(yNew,i);
        else
            i = i + 1;
        end
    else
        i = i + 1;
    end
    
    iMax = length(sNew);
    
end

    %Insert a new column at the desired index
    function newArray = ColumnInsert(val,arrayIn,index)
        A = arrayIn(:,1:(index-1));
        B = arrayIn(:,index:end);
        newArray = [A,val,B];
    end

    %Remove the column at the desired index
    function newArray = ColumnRemove(arrayIn,index)
        A = arrayIn(:,1:(index-1));
        B = arrayIn(:,(index+1):end);
        newArray = [A,B];
    end

end