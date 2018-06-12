function [points state] = model_sampling_ACHR(model,warmupPoints, pointsCount, stepsPerPoint, state)
% ACHRSampler Artificial Centering Hit-and-Run sampler
%
% ACHRSampler(model,warmupPoints,fileName,nFiles,pointsPerFile,stepsPerPoint,initPoint)
%
% model    Model structure
% warmupPoints Warmup points
% pointsCount Number of sample points
% stepsPerPoint Number of sampler steps per point saved
% initPoint     Initial point (optional)
%
% Markus Herrgard, Gregory Hannum, Ines Thiele, Nathan Price 4/14/06

warning off MATLAB:divideByZero;

if isfield(model, 'N')
    N = model.N;
else
    N = null(full(model.S));
end

% Minimum allowed distance to the closest constraint
maxMinTol = 1e-9;
% Ignore directions where u is really small
uTol = 1e-9; 
% Project out of directions that are too close to the boundary
dTol = 1e-14;

% Number of warmup points
[nRxns, nWrmup] = size(warmupPoints);


% Set the start point
if (nargin < 5) || isempty(state)
    % Find the center of the space
    centerPoint = mean(warmupPoints,2);
    prevPoint = centerPoint;
    totalStepCount = 0;
else
    centerPoint = state.centerPoint;
    prevPoint = state.prevPoint;
    totalStepCount = state.totalStepCount;
end


fprintf('start sampling\n');
% Allocate memory for all points
points = zeros(nRxns,pointsCount); 

pointCount = 1;
while pointCount <= pointsCount

    % Create the random step size vector
    randVector = rand(stepsPerPoint,1);

    stepCount = 1;
    while (stepCount <= stepsPerPoint)

        % Pick a random warmup point
        randPointID = ceil(nWrmup*rand);
        randPoint = warmupPoints(:,randPointID);

        % Get a direction from the center point to the warmup point
        u = randPoint - centerPoint;
        u = u / norm(u);

        % Figure out the distances to upper and lower bounds
        distUb = model.ub - prevPoint;
        distLb = prevPoint - model.lb;

        % Figure out if we are too close to a boundary
        validDir = (distUb > dTol) & (distLb > dTol);

        % Figure out positive and negative directions
        validU = u(validDir);
        posDirn = validU > uTol;
        negDirn = validU < -uTol;

        % Figure out all the possible maximum and minimum step sizes
        inverseValidU = 1 ./ validU;
        maxStepTemp = distUb(validDir) .* inverseValidU;
        minStepTemp = -distLb(validDir) .* inverseValidU;
        maxStepVec = [maxStepTemp(posDirn); minStepTemp(negDirn)];
        minStepVec = [minStepTemp(posDirn); maxStepTemp(negDirn)];

        % Figure out the true max & min step sizes
        maxStep = min(maxStepVec);
        minStep = max(minStepVec);
        %fprintf('%f\t%f\n',minStep,maxStep);

        % Find new direction if we're getting too close to a constraint
        if (abs(minStep) < maxMinTol && abs(maxStep) < maxMinTol) || (minStep > maxStep)
            fprintf('Warning %f %f\n',minStep,maxStep);
            continue;
        end

        % Pick a rand out of list_of_rands and use it to get a random
        % step distance
        stepDist = randVector(stepCount)*(maxStep-minStep)+minStep;

        % Advance to the next point
        curPoint = prevPoint + stepDist*u;

        % Reproject the current point and go to the next step
        if mod(totalStepCount,10) == 0
            if (full(max(max(abs(model.S*curPoint)))) > 1e-9)
              curPoint = N*(N'*curPoint);
            end
        end

        overInd = model.ub < curPoint;
        underInd = model.lb > curPoint;
        if any(overInd) || any(underInd)
            curPoint(overInd) = model.ub(overInd);
            curPoint(underInd) = model.lb(underInd);
        end


        prevPoint = curPoint;
        stepCount = stepCount + 1;

        % Count the total number of steps
        totalStepCount = totalStepCount + 1;
        if mod(totalStepCount, 1000) == 0
            print_progress((stepCount+pointCount*stepsPerPoint) / (pointsCount * stepsPerPoint));
        end
        %recalculate the center point
        centerPoint = ((nWrmup+totalStepCount)*centerPoint + curPoint)/(nWrmup+totalStepCount+1);

    end % Steps per point

    % Add the current point to points
    points(:,pointCount) = curPoint;

    pointCount = pointCount + 1;

end % Points
state.centerPoint = centerPoint;
state.prevPoint = prevPoint;
state.totalStepCount = totalStepCount;
fprintf('done sampling\n');
