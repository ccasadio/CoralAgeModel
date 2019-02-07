numPeriods = 25;
pointsPerPeriod = 12;
upscaleMult = 25;
samplingRate = 1.2; % 1.2cm/yr
numRealizations = 1000;

%% construct a growth rate series 
lagCoef1 = 0.2378;
lagCoef2 = 0.2022;
growthRateSigma = .2;
grModel = arima('Constant',0,'AR',{lagCoef1, lagCoef2},'Variance',growthRateSigma*growthRateSigma);
% evaluate the AR model and make the mean 1.2 cm per year
[rateOfGrowth,~] = simulate(grModel, numPeriods);
rateOfGrowth = rateOfGrowth + (samplingRate - mean(rateOfGrowth));

% change to monotonically increasing depth vector
rateOfGrowth = cumsum([0,rateOfGrowth']);

% interpolate the AR model to increase the sample rate by the upscaleMult
rateOfGrowth = interp1(1:length(rateOfGrowth), rateOfGrowth, linspace(1, length(rateOfGrowth), length(rateOfGrowth)*upscaleMult*pointsPerPeriod), 'linear');

%% create signal and add gaussian noise (signal "fingerprint")
signal = sin(2*pi*linspace(0,numPeriods,length(rateOfGrowth)));
signal = signal + randn(size(signal))*.1;

%% average signal (sampling model)
% ideal agemodel result
for i = 1:(numPeriods*pointsPerPeriod + 1)
    signal_ideal(i) = mean(signal( ((i-1)*upscaleMult+1):(i*upscaleMult) ));
end

% sampled signal
% use rate of growth and constant sampling rate
for i = 1:(numPeriods*pointsPerPeriod + 1)
    signal_sampled(i) = mean(signal( (rateOfGrowth > ((i-1)*samplingRate/pointsPerPeriod)) & (rateOfGrowth < (i*samplingRate/pointsPerPeriod)) ));
end

%% montecarlo agemodel to see if noise can help converge on a better solution and which correctness estimates are most useful
x = 0:(1/pointsPerPeriod):numPeriods; 
sensitivity = .7:.01:1.2;
results = cell(numRealizations,length(sensitivity));
correctness_est = cell(numRealizations,length(sensitivity));
correctness = cell(numRealizations,length(sensitivity));
valSen = cell(numRealizations,length(sensitivity));
for r = 1:numRealizations
    % add gaussian noise
    y = signal_sampled + randn(size(signal_sampled))*.1;

    %% iterate through spline sensitivities
    for s = 1:length(sensitivity)
        valSen{r,s} = sensitivity(s);
        [ts, criticalPoints] = corallinagemodel([x;y]','sensitivity',sensitivity(s));
        results{r,s} = ts(:,2)';
            
        %% calculate estimated correctness metrics

        %% calculate actual correctness
        l = min([length(results{r,s}),length(signal_ideal)]);
        lDiff = abs(length(results{r,s})-length(signal_ideal));
        err = results{r,s}(1:l) - signal_ideal(1:l);
        correctness{r,s} = mean(([err,max(err)*ones(1,lDiff)]).^2 );
    end
end

%% plot results
figure
boxplot([correctness{:}],[valSen{:}])

for r = 1:numRealizations
    for s = 1:length(sensitivity)
        valLen{r,s} = length(signal_ideal)-length(results{r,s});
    end
end


figure
boxplot([valLen{:}],[valSen{:}])
