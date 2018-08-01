function [ts, criticalPoints] = corallinagemodel( inputData, varargin )
%CORALLINAGEMODEL Summary of this function goes here
%   Detailed explanation goes here

defaultSensitivity = 1;
p = inputParser;

addRequired(p, 'inputData', @(x) isnumeric(x) && ismatrix(x) && ((size(x,2) == 2) || size(x,1) == 2 ));
addParameter(p, 'sensitivity', defaultSensitivity, @isnumeric);
addParameter(p, 'numyears', 0, @isnumeric);
addParameter(p, 'startperiodicity', 12, @isnumeric);
addParameter(p, 'endperiodicity', 12, @isnumeric);

parse(p,inputData,varargin{:});

if size(inputData, 1) == 2
   inputData = inputData'; 
end

% find and remove NaNs
inputData = inputData( ~(isnan(inputData(:,1)) | isnan(inputData(:,2))),:);

[pow,w] = pwelch(inputData(:,2));
[~, loc] = findpeaks(pow, w, 'SortStr', 'descend');

% convert to period
loc = (2*pi)./loc;

try
    loc = loc(loc < (p.Results.startperiodicity * 1.5));
    loc = loc(loc > p.Results.startperiodicity * 0.5);
    ppy = loc(1);
catch
    %default to 12 ppy if power spectrum fails to pick the correct amount
    ppy = p.Results.startperiodicity;
end

[ts, criticalPoints] = linearAgeModel(inputData, ppy, p.Results.sensitivity*2.7, p.Results.endperiodicity);

% if a year target is set then optimize sensitivity to get that number of
% years
if p.Results.numyears > 0
    
    i = 0;
    delta = .1;
    proximity = p.Results.numyears - ts(end,1);
    sensitivity = p.Results.sensitivity;
    
    tolerance = 1;
    maxIterations = 25;
    bounds = [-inf inf];

    %binary search algorith
    while (i < maxIterations) && (abs(proximity) > tolerance)
        [ts, criticalPoints] = linearAgeModel(inputData, ppy, sensitivity*2.7, p.Results.endperiodicity);
        proximity = p.Results.numyears - ts(end,1);

        if proximity > 0
            if bounds(1) < sensitivity
                bounds(1) = sensitivity;
            end
        else
            if bounds(2) > sensitivity
                bounds(2) = sensitivity;
            end
        end

        if (sum(abs(bounds)) == inf)

            if sum(bounds) > 0
                sensitivity = sensitivity + delta;
            else
                sensitivity = sensitivity - delta;
            end

            delta = delta*2;

        else
            sensitivity = mean(bounds);
        end
        i = i + 1;
    end
    
    if i == 25
        warning('Max number of iterations (%d) reached.',i);
    end
end

end


function [ts, criticalPoints] = linearAgeModel(data, pointsPerYear, splineSensitivity, requiredPointsPerYear)

debug = false;

% Create new interpolated timeseries
%  Constructing a cubic spline fit to the data using numknots (# of knots)
%  = the number of bins * an arbitrary splinesensitivity input above
numknots = (length(data)/pointsPerYear)*splineSensitivity;

warning('off','SPLINES:CHCKXYWP:NaNs');
sp = spap2(numknots,3,data(:,1),data(:,2));

if debug
    figure;
    plot(data(:,1),fnval(sp,data(:,1)),'k-');
    hold on;
    plot(data(:,1),data(:,2),'ro');
    hold off;
end

% Calculating the first derivative of the spline fit
firstder = fnder(sp,1);

% inflection points of the spline fit to the original data
criticalPointsSpline = fnzeros(firstder);
criticalPointsSpline = mean(criticalPointsSpline,1)';

% preallocate
intermediateX = zeros((length(criticalPointsSpline)-1)*6, 1);

% create x values that are evenly spaced between critical points
for i = 1:(length(criticalPointsSpline)-1)
    intermediateX((1+((i-1)*6)):(1+(i*6)),1) = linspace(criticalPointsSpline(i,1),criticalPointsSpline(i+1,1), requiredPointsPerYear/2 + 1);
end

% add x values for endpoints
meanLength = mean(diff(intermediateX(:,1)));
intermediateX = [(data(1,1):meanLength:intermediateX(1,1))'; intermediateX(:,1); (intermediateX(end,1):meanLength:data(end,1))'];

% linearly interpolate data to 12 ppy
ts(:,2) = interp1(data(:,1), data(:,2), intermediateX);
ts(:,1) = 0:(1/requiredPointsPerYear):(length(ts(:,2))-1)/requiredPointsPerYear;
criticalPoints = criticalPointsSpline;
end
