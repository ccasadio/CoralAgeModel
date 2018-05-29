function [out] = corallinagemodel( inputData )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% find and remove NaNs
inputData = inputData( ~(isnan(inputData(:,1)) | isnan(inputData(:,2))),:);

[pow,w] = pwelch(inputData(:,2));
[~, loc] = findpeaks(pow, w, 'SortStr', 'descend');

% convert to frequency
loc = (2*pi)./loc;

try
    loc = loc(loc < 18);
    loc = loc(loc > 6);
    ppy = loc(1);
catch
    %default to 12 ppy if power spectrum fails to pick the correct amount
    ppy = 12;
end

out = linearAgeModel(inputData,ppy,3.7);
end

function [ts, criticalPoints] = linearAgeModel(data, pointsPerYear, splineSensitivity)
requiredPointsPerYear = 12;

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
