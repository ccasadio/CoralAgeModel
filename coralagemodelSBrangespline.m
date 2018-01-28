

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USER INPUTS  %

% User enters the number of analyses per year
% function of growth-rate and sampling interval
numberpointsperyear = 12; %15 for TMI1, 10 for TMS5

% Another user input for the sensitivity of the spline fit. Closer to 2
% means a lower fit. Closer to 3 means a better fit - but caution that the
% closer to 3 may result in spurious years.
splinesensitivity = 2.6; %2.6 for TMI1, 2.7 for TMS5

x = (1:3.1/12:110)';
y = (sin(x) + rand(size(x))./4)./7 + 7.7;

ts = ModCoral_06SBA1_Depth_SrCa;

num = numberpointsperyear - 1;

numbins = round(length(ts)/numberpointsperyear);

% create new interpolated timeseries
% Constructing a cubic spline fit to the data using numknots (# of knots)
% = the number of bins * an arbitrary splinesensitivity input above
numknots = numbins*splinesensitivity;

drilldepth = ts(:,1);
d18Odepth = ts(:,2);
weight = double(~(isnan(d18Odepth)));

sp = spap2(numknots,3,drilldepth,d18Odepth,weight);
splined18O = fnval(sp,drilldepth);


% Calculating the first derivative of the spline fit
firstder = fnder(sp,1);
fderval = fnval(firstder,drilldepth);

% Figure to plot the spline fit, original timeseries and the first
% derivative to determine how well the inflection points look to be near
% the summers and winters
figure(1)
subplot(2,1,1)
hold on
plot(drilldepth,fderval)
line([0 length(drilldepth)],[0 0],'Color','k')
set(gca,'Xlim',[0 max(ts(:,1)+10)])
subplot(2,1,2)
hold on
set(gca,'YDir','reverse','Xlim',[0 max(ts(:,1)+10)])
plot(drilldepth,splined18O,'b-')
plot(ts(:,1),ts(:,2),'o')

% Logging the inflection points where the first derivative = 0 (sign
% change)
counter = 1;
for j = 1:length(fderval)-1
    if fderval(j,1) > 0 && fderval(j+1,1) < 0
        if isnan(ts(j,2)) == 1
            inflections(counter,:) = [ts(j,1) splined18O(j)];
        else
            inflections(counter,:) = ts(j,:);
        end
        counter = counter + 1;
    elseif fderval(j,1) < 0 && fderval(j+1,1) > 0
        if isnan(ts(j,2)) == 1
            inflections(counter,:) = [ts(j,1) splined18O(j)];
        else
            inflections(counter,:) = ts(j,:);
        end
        counter = counter + 1;
    else
        counter = counter;
    end
end

% computing mean "distance" for a year. This number should be close to the
% number input by the user in the beginning of the program
for ii = 1:length(inflections)-1
    meanyears(ii,1) = inflections(ii+1,1)-inflections(ii,1);
end
meandistance = mean(meanyears);

% Number of years in the record which is half of the number of inflections
% points
numberofyearsinrecord = length(inflections)/2;

% computing seasonal cycle
% compiling all of the winters, which is the max value in the
% surrounding five data points from the inflection point
% If data gaps are filled using NaN, then the splined18O value is used
counterer = 1;
for j = 2:length(fderval)-1
    if fderval(j,1) > 0 && fderval(j+1,1) < 0
        if isnan(ts(j-2,2))||isnan(ts(j-1,2))||isnan(ts(j,2))||isnan(ts(j+1,2))||isnan(ts(j+2,2)) == 1;
            if isnan(ts(j-2,2)) == 1
                newminus2 = [ts(j-2,1) splined18O(j-2,:)];
            else
                newminus2 = ts(j-2,:);
            end
            if isnan(ts(j-1,2)) == 1
                newminus1 = [ts(j-1,1) splined18O(j-1,:)];
            else
                newminus1 = ts(j-1,:);
            end
            if isnan(ts(j,2)) == 1
                newone = [ts(j,1) splined18O(j,:)];
            else
                newone = ts(j,:);
            end
            if isnan(ts(j+1,2)) == 1
                newplus1 = [ts(j+1,1) splined18O(j+1,:)];
            else
                newplus1 = ts(j+1,:);
            end
            if isnan(ts(j+1,2)) == 1
                newplus2 = [ts(j+2,1) splined18O(j+2,:)];
            else
                newplus2 = ts(j+2,:);
            end
            winterblock = [newminus2;newminus1;newone;newplus1;newplus2];
            [r,c,v] = find(winterblock == max(winterblock(:,2)));
            winter(counterer,:) = winterblock(r(1,1),:);
        else
            winterblock = [ts(j-2,:);ts(j-1,:);ts(j,:);ts(j+1,:);ts(j+2,:)];
            [r,c,v] = find(winterblock == max(winterblock(:,2)));
            winter(counterer,:) = winterblock(r(1,1),:);
        end
        counterer = counterer + 1;
    else
        counterer = counterer;
    end
end

% compiling all of summers, which is the min value in the
% surrounding five data points from the inflection point
% If data gaps are filled using NaN, then the splined18O value is used
counterer = 1;
for j = 2:length(fderval)-1
    if fderval(j,1) < 0 && fderval(j+1,1) > 0
        if isnan(ts(j-2,2))||isnan(ts(j-1,2))||isnan(ts(j,2))||isnan(ts(j+1,2))||isnan(ts(j+2,2)) == 1;
            if isnan(ts(j-1,2)) == 1
                newminus2 = [ts(j-2,1) splined18O(j-2,:)];
            else
                newminus2 = ts(j-2,:);
            end
            if isnan(ts(j-1,2)) == 1
                newminus1 = [ts(j-1,1) splined18O(j-1,:)];
            else
                newminus1 = ts(j-1,:);
            end
            if isnan(ts(j,2)) == 1
                newone = [ts(j,1) splined18O(j,:)];
            else
                newone = ts(j,:);
            end
            if isnan(ts(j+1,2)) == 1
                newplus1 = [ts(j+1,1) splined18O(j+1,:)];
            else
                newplus1 = ts(j+1,:);
            end
            if isnan(ts(j+2,2)) == 1
                newplus2 = [ts(j+2,1) splined18O(j+2,:)];
            else
                newplus2 = ts(j+2,:);
            end
            summerblock = [newminus2;newminus1;newone;newplus1;newplus2];
            [r,c,v] = find(summerblock == min(summerblock(:,2)));
            summer(counterer,:) = summerblock(r(1,1),:);
        else
            summerblock = [ts(j-2,:);ts(j-1,:);ts(j,:);ts(j+1,:);ts(j+2,:)];
            [r,c,v] = find(summerblock == min(summerblock(:,2)));
            summer(counterer,:) = summerblock(r(1,1),:);
        end
        counterer = counterer + 1;
    else
        counterer = counterer;
    end
end

% If one manually changes the summers and winters then run from below to
% the end to recalculate the seasonal cycle calculations
%%
% Uncomment the lines below if the summers and winters are manually changed
% then run just this cell to recalculate the seasonal cycle calc

% close all
% clear  newdepth scnolag plotnolag sclag plotlag coralinterp coralts time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolating the timeseries to equal monthly spacing using the peaks and
% valleys, or summers and winters picked previously to make 5 datapoints
% between each peak and valley to have 12 points/months in a year

% Initialize timeseries for spline (depth, time, var)
splinets = zeros(length(summer) + length(winter),3);

counttime = 1;
othercounttime = 1;
% Case 1
if summer(1,1) < winter (1,1) && length(summer) > length(winter)
    for i = 1:length(summer)-1
        splinets(counttime,1) = summer(i,1);
        splinets(counttime+1,1) = winter(i,1);
        splinets(counttime,2) = othercounttime-1+0.625;
        splinets(counttime,3) = summer(i,2);
        splinets(counttime+1,3) = winter(i,2);
        splinets(counttime+1,2) = othercounttime-1+1.125;
        counttime = counttime + 2;
        othercounttime = othercounttime + 1;
    end
    for i = length(summer)-1:length(summer)
        splinets(counttime,1) = summer(i,1);
        splinets(counttime,2) = othercounttime-1+0.625;
        splinets(counttime,3) = summer(i,2);
    end
    % Case 2
elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
    for i = 1:length(summer)
        splinets(counttime,1) = summer(i,1);
        splinets(counttime+1,1) = winter(i,1);
        splinets(counttime,2) = othercounttime-1+0.625;
        splinets(counttime,3) = summer(i,2);
        splinets(counttime+1,3) = winter(i,2);
        splinets(counttime+1,2) = othercounttime-1+1.125;
        counttime = counttime + 2;
        othercounttime = othercounttime + 1;
    end
    % Case 3
elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
    for i = 1:length(winter)-1
        splinets(counttime,1) = winter(i,1);
        splinets(counttime+1,1) = summer(i,1);
        splinets(counttime,2) = othercounttime-1+0.125;
        splinets(counttime,3) = winter(i,2);
        splinets(counttime+1,3) = summer(i,2);
        splinets(counttime+1,2) = othercounttime-1+0.625;
        counttime = counttime + 2;
        othercounttime = othercounttime + 1;
    end
    for i = length(winter)-1:length(winter)
        splinets(counttime,1) = winter(i,1);
        splinets(counttime,2) = othercounttime-1+0.625;
        splinets(counttime,3) = winter(i,2);
    end
    % Case 4
elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
    for i = 1:length(winter)-1
        splinets(counttime,1) = winter(i,1);
        splinets(counttime+1,1) = summer(i,1);
        splinets(counttime,2) = othercounttime-1+0.125;
        splinets(counttime,3) = winter(i,2);
        splinets(counttime+1,3) = summer(i,2);
        splinets(counttime+1,2) = othercounttime-1+0.625;
        counttime = counttime + 2;
        othercounttime = othercounttime + 1;
    end
end


coralts(:,1) = [min(splinets(:,2)):1/12:max(splinets(:,2))]';
coralts(:,2) = interp1(splinets(:,2),splinets(:,3),coralts(:,1),'pchip');

% timestarts on the month of the first seasonal extreme and ends on the
% month of the last seasonal extreme

time = coralts(:,1);

% scnolag = the first summer or winter is difference from the next;
% if the seasons are not the same length, then...

ct = 1;
% Case 1
if summer(1,1) < winter (1,1) && length(summer) > length(winter)
    scnolag=zeros(length(winter(:,1)),1);
    plotnolag=zeros(length(winter(:,1)),1);
    sclag=zeros(length(winter(:,1)),1);
    plotlag=zeros(length(winter(:,1)),1);
    for k = 1:2:length(winter)*2-1
        scnolag(ct,1) = splinets(k,3)-splinets(k+1,3);
        plotnolag(ct,1) = (splinets(k+1,1)+splinets(k,1))/2.;
        sclag(ct,1) = splinets(k+2,3)-splinets(k+1,3);
        plotlag(ct,1) = (splinets(k+2,1)+splinets(k+1,1))/2.;
        ct = ct + 1;
    end
    % Case 2
elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
    scnolag=zeros(length(winter(:,1)),1);
    plotnolag=zeros(length(winter(:,1)),1);
    sclag=zeros(length(winter(:,1)),1);
    plotlag=zeros(length(winter(:,1)),1);
    for k = 1:2:length(winter)*2-1
        scnolag(ct,1) = splinets(k,3)-splinets(k+1,3);
        plotnolag(ct,1) = (splinets(k+1,1)+splinets(k,1))/2.;
        sclag(ct,1) = splinets(k+2,3)-splinets(k+1,3);
        plotlag(ct,1) = (splinets(k+2,1)+splinets(k+1,1))/2.;
        ct = ct + 1;
    end
    % Case 3
elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
    scnolag=zeros(length(summer(:,1)),1);
    plotnolag=zeros(length(summer(:,1)),1);
    sclag=zeros(length(summer(:,1)),1);
    plotlag=zeros(length(summer(:,1)),1);
    for k = 1:2:length(summer)*2-1
        scnolag(ct,1) = splinets(k+1,3)-splinets(k,3);
        plotnolag(ct,1) = (splinets(k+1,1)+splinets(k,1))/2.;
        sclag(ct,1) = splinets(k+2,3)-splinets(k+1,3);
        plotlag(ct,1) = (splinets(k+2,1)+splinets(k+1,1))/2.;
        ct = ct + 1;
    end
    % Case 4
elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
    scnolag=zeros(length(summer(:,1)),1);
    plotnolag=zeros(length(summer(:,1)),1);
    sclag=zeros(length(summer(:,1)),1);
    plotlag=zeros(length(summer(:,1)),1);
    for k = 1:2:length(summer)*2-1
        scnolag(ct,1) = splinets(k+1,3)-splinets(k,3);
        plotnolag(ct,1) = (splinets(k+1,1)+splinets(k,1))/2.;
        sclag(ct,1) = splinets(k+2,3)-splinets(k+1,3);
        plotlag(ct,1) = (splinets(k+2,1)+splinets(k+1,1))/2.;
        ct = ct + 1;
    end
end

ssncycnolag = mean(scnolag);
ssncyclag = mean(sclag);

% Figure to check if the first derivative picks summers and winters
% correctly where the original and spline timeseries are also plotted
figure(2)
hold on
set(gca,'YDir','reverse','Xlim',[0 max(ts(:,1)+10)])
plot(drilldepth,splined18O,'b-')
plot(ts(:,1),ts(:,2),'k-')
plot(summer(:,1),summer(:,2),'ro')
plot(winter(:,1),winter(:,2),'bo')
for i = 1:length(summer)
    plotsummer = summer(i);
    plots = num2str(plotsummer);
    text(summer(i,1),summer(i,2)-0.1,(plots),'Color','r');
end
for i = 1:length(winter)
    plotwinter = winter(i);
    plotw = num2str(plotwinter);
    text(winter(i,1),winter(i,2)+0.1,(plotw),'Color','b');
end
hold off

figure(3)
subplot(2,1,1)
hold on
set(gca,'YDir','reverse','Xlim',[0 max(ts(:,1)+10)])
plot(ts(:,1),ts(:,2),'k-')
hold off
subplot(2,1,2)
hold on
set(gca,'YDir','reverse','Xlim',[0 max(time)+1])
plot(coralts(:,1),coralts(:,2))
hold off
%%
figure(4)
subplot(3,1,1)
hold on
plot((2008-numberofyearsinrecord+(drilldepth/12)),splined18O,'b-')
plot((2008-numberofyearsinrecord+(ts(:,1)/12)),ts(:,2),'k-')
plot((2008-numberofyearsinrecord+(summer(:,1)/12)),summer(:,2),'ro')
plot((2008-numberofyearsinrecord+(winter(:,1)/12)),winter(:,2),'bo')
set(gca,'YDir','reverse');%,'Xlim',[1845 2008],'xtick',[1850:10:2000])
xlabel('Year')
ylabel('Coral d18O (permil)')
hold off
subplot(3,1,2)
hold on
bar(plotnolag,-1*scnolag,1,'k')
line([0 max(ts(:,1)+10)],[-1*ssncycnolag -1*ssncycnolag])
set(gca,'Xlim',[-30 max(ts(:,1)+10)])
ylabel('Seasonal Cycle Amplitude')
hold off
subplot(3,1,3)
hold on
bar(plotlag,-1*sclag,1,'g')
line([0 max(ts(:,1)+10)],[-1*ssncyclag -1*ssncyclag])
set(gca,'Xlim',[-30 max(ts(:,1)+10)])
ylabel('Seasonal Cycle Amplitude')
hold off

% Doing a bootstrap analysis to determine the 95% confidence intervals for
% the mean of the computed climatologies
scmeannolag = @(scnolag)mean(scnolag);
scmeanlag = @(sclag)mean(sclag);

nolagbootci = bootci(1000,scmeannolag,scnolag);
lagbootci = bootci(1000,scmeanlag,sclag);

avgsc = -1*(ssncycnolag + ssncyclag)/2.;
avgerror = -1*((nolagbootci(1,1) - ssncycnolag) + (ssncycnolag - nolagbootci(2,1)) + (lagbootci(1,1) - ssncyclag) + (ssncyclag - lagbootci(2,1)))/4.;


