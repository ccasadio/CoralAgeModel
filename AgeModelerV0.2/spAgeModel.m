function [ts, criticalPoints] = spAgeModel(data, pointsPerYear, splineSensitivity)
requiredPointsPerYear = 12;

% Create new interpolated timeseries
%  Constructing a cubic spline fit to the data using numknots (# of knots)
%  = the number of bins * an arbitrary splinesensitivity input above
numknots = (length(data)/pointsPerYear)*splineSensitivity;

warning('off','SPLINES:CHCKXYWP:NaNs');
sp = spap2(numknots,3,data(:,1),data(:,2));

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