function ylist = MyTrackCurve(userf,y0,ytan,varargin)
%MyTrackCurve Pseudo-Arclength continuation of a family of input functions with variable stepsize
%   Input:
%   
%   userf - Column vector, function to track
%   y0 - Column vector, numeric, starting point of arc
%   ytan - Estimate for initial tangent
%   varargin - Optional arguments
%
%   Output:
%   
%   ylist - tracked points
%

parser = inputParser;

addRequired(parser,'userf',@(func) isa(func,'function_handle'));
addRequired(parser,'y0',@isnumeric);
addRequired(parser,'ytan',@isnumeric);
addOptional(parser,'stepsize',1e-2,@isnumeric);                                 %Stepsize of first guess
addOptional(parser,'maxSolveIter',10,@isnumeric);                               %Maximum number of MySolve iterations
addOptional(parser,'nMax',100,@isnumeric);                                      %Maximum number of iterations for curve tracking
addOptional(parser,'jStep',1e-8,@isnumeric);                                    %Jacobian accuracy
addOptional(parser,'stop',@(y) false,@(func) isa(func,'function_handle'));      %Condition to stop curve tracking
addOptional(parser,'sMax',1e-1,@isnumeric);                                     %Maximum allowed step size
addOptional(parser,'sMin',1e-12,@isnumeric);                                    %Minimum allowed step size
addOptional(parser,'doPrint',false,@islogical);                                 %Print certain alerts or warnings
addOptional(parser,'correctGuess',false,@islogical);                            %Attempt to converge user's first guess to a true root
addOptional(parser,'tol',1e-6,@isnumeric);                                      %MySolve tolerance
parse(parser,userf,y0,ytan,varargin{:})

stepsize = parser.Results.stepsize;
maxSolveIter = parser.Results.maxSolveIter;
nMax = parser.Results.nMax;
jStep = parser.Results.jStep;
stop = parser.Results.stop;
sMax = parser.Results.sMax;
sMin = parser.Results.sMin;
doPrint = parser.Results.doPrint;
correctGuess = parser.Results.correctGuess;
solveTol = parser.Results.tol;

%Intialise output matrix
ylist = NaN(length(y0),nMax);

i = 1;
yk = y0;

while i <= nMax
        
    %If correctGuess is true, converge first guess to closes true root
    if (i == 1 && correctGuess)
        stepsize = 0;
    elseif (i == 2 && correctGuess)
        stepsize = parser.Results.stepsize;
    end

    %Step forward in tangent direction
    yp = yk + stepsize.*ytan;
    
    %Set new functions for arclength continuation
    f = @(y) [userf(y); ytan.' * (y-yp)];
    df = @(y) MyJacobian(f,y,jStep);
       
    %Solve at new point
    [ykn,converged,~] = MySolve(f,yp,df,'tol',solveTol,'maxIter',maxSolveIter);
        
    while(converged == 0)
        
        %If convergence failed, try a smaller step, repeating above proceedure
        stepsize = stepsize/2;
        
        if(stepsize < sMin) %Convergence fails if stepsize is below lower limit
            warning("Failed to converge - Step size too small. This issue may potentially be solved by setting correctGuess to false.")
            return
        end
        
        yp = yk + stepsize.*ytan;
        
        f = @(y) [userf(y); ytan.' * (y-yp)];
        df = @(y) MyJacobian(f,y,jStep);
        
        [ykn,converged,~] = MySolve(f,yp,df,'tol',solveTol,'maxIter',maxSolveIter);
        
    end
    
    %Accept new root
    yk = ykn;
      
    %Store new root
    ylist(:,i) = yk;
    
    %Test stop condition, and break if met
    if(stop(yk))
        break
    end     
    
    %Calculate new tangent
    A = [MyJacobian(userf,yk,jStep);ytan.'];
    I = eye(length(A));
    B = I(:,end);
    
    z=A\B;
    ytan = (z./norm(z,Inf)) * sign(z' * ytan);
        
    %Increase stepsize for next iteration
    stepsize = min(stepsize*1.25,sMax);
    
    i = i + 1;
    
end

%Suggest improved nMax value to user
if (i < nMax && doPrint)
    disp(strcat("MyTrackCurve converged after ", num2str(i), " iterations. Consider using a smaller value of nMax to save memory."))
end

end

