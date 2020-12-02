function [yNew,xNew,sNew] = MapLine(M,x,s,maxDist,maxAngle)

yNew = M(x);
sNew = s;
xNew = x;

iMax = length(sNew);
i = 1;
while i < iMax
    
    while( norm(yNew(:,i+1) - yNew(:,i)) > maxDist)
        
        sStar = (sNew(i) + sNew(i+1)) ./ 2;
        xStar = interp1(sNew,xNew',sStar,'spline')';
        yStar = M(xStar);
        
        sNew = ColumnInsert(sStar,sNew,i+1);
        xNew = ColumnInsert(xStar,xNew,i+1);
        yNew = ColumnInsert(yStar,yNew,i+1);
        
    end
    
    if(i ~= 1)
               
        segA = yNew(:,i) - yNew(:,i-1);
        segB = yNew(:,i+1) - yNew(:,i-1);
        segAngle = acos(dot(segA./norm(segA),segB./norm(segB)));
        
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
    
    segA = yNew(:,i) - yNew(:,i-1);
    segB = yNew(:,i+1) - yNew(:,i-1);
    segAngle = acos(dot(segA./norm(segA),segB./norm(segB)));
   
    if(segAngle < maxAngle)            
        if(norm(yNew(:,i+1) - yNew(:,i-1)) < maxDist)
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

    function newArray = ColumnInsert(val,arrayIn,index)
        A = arrayIn(:,1:(index-1));
        B = arrayIn(:,index:end);
        newArray = [A,val,B];
    end

    function newArray = ColumnRemove(arrayIn,index)
        A = arrayIn(:,1:(index-1));
        B = arrayIn(:,(index+1):end);
        newArray = [A,B];
    end

end