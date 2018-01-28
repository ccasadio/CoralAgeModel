clear all
close all
warning off
format long

%Dock all figures by default
set(0,'DefaultFigureWindowStyle','docked');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USER INPUTS  %

% choose either 2 for d18O, 3 for Sr/Ca or 4 for d13C
variable = 3;

% User enters the number of analyses per year
% function of growth-rate and sampling interval
numberpointsperyear = 15; %10 points per year for TMS5; 15 for TMI1

% Another user input for the sensitivity of the spline fit. Closer to 2
% means a lower fit. Closer to 3 means a better fit - but caution that the
% closer to 3 may result in spurious years.
splinesensitivity = 2.6; %2.7 for TMS5, %2.6 for TMI1

% Reading in data from the master sheet where
% column 1 = sample #
% column 2 = d18O
% column 3 = Sr/Ca
% column 4 = d13C
% Uncomment whichever coral you want to run. Uncommenting them all makes
% the program run slow.
%
%open data for fossil corals 11-TM-S5 and 11-TM-I1
filename = '/Users/alawman/Documents/TM Last 2k Coral Analysis/VanuatuCoralDataComposites.xlsx';
%sheet = '11-TM-S5';
sheet = '11-TM-I1';
%range = 'C4:F1194'; %TM-S5
range = 'C249:F576'; %TM-I1 (Paths A-C only)
[numCoral,txt] = xlsread(filename, sheet, range);

 TMI1 = [];
 TMI1 = numCoral;
% TMI1 = xlsread('TMIFossilCoralInput.xlsx');
% TMI1 = xlsread('TMIFossilCoralInput.xlsx','TMI1_depthfrombottom');
% TMS5 = []; %use 2.7 spline sensitivity and 10 points per year for TMS5
% TMS5 = num_TMS5;
% TMS5 = xlsread('TMS5FossilCoralInput.xlsx','11TMS5_2');
% TMS5 = xlsread('TMS5FossilCoralInput.xlsx','11TMS5_3');
% TMS5 = xlsread('TMS5_FossilCoralInput_suboptimalComposite.xlsx','TMS5_suboptimal');
% TMN = [];
% TMN = xlsread('TMNFossilCoralInput.xlsx','05TMN');
% SB = [];
% SB = xlsread('SBFossilCoralInput.xlsx','SBrawSrCa');
% r98k3 = [];
% r98k3 = xlsread('FossilCoralInput.xls','r98k3');
% SAG = [];
% SAG = xlsread('FossilCoralInput.xls','SAG');
% Junhua = [];
% Junhua = xlsread('FossilCoralInput.xls','Junhua');
% T01_B_Mo = [];
% T01_B_Mo = xlsread('FossilCoralInput.xls','T01_B_Mo');
% T01_B =[];
% T01_B = xlsread('FossilCoralInput.xls','T01_B');
% RND99_Mo = [];
% RND99_Mo = xlsread('FossilCoralInput.xls','RND99_Mo');
% RND99 = [];
% RND99 = xlsread('FossilCoralInput.xls','RND99');
% T01_AJ_4 = [];
% T01_AJ_4 = xlsread('FossilCoralInput.xls','T01_AJ_4');
% T01_AD_1 = [];
% T01_AD_1 = xlsread('FossilCoralInput.xls','T01_AD_1');
% T01_Z_1a = [];
% T01_Z_1a = xlsread('FossilCoralInput.xls','T01_Z_1a');
% T01_AP_1_2 = [];
% T01_AP_1_2 = xlsread('FossilCoralInput.xls','T01_AP_1_2');
% T01_AQ_Depth = [];
% T01_AQ_Depth = xlsread('FossilCoralInput.xls','T01_AQ_Depth');
% T01_AQ_Mo = [];
% T01_AQ_Mo = xlsread('FossilCoralInput.xls','T01_AQ_Mo');
% T01_AR_1_2 = [];
% T01_AR_1_2 = xlsread('FossilCoralInput.xls','T01_AR_1_2');
% T01_D_1_2 = [];
% T01_D_1_2 = xlsread('FossilCoralInput.xls','T01_D_1_2');
% TET99_B_Mo = [];
% TET99_B_Mo = xlsread('FossilCoralInput.xls','TET99_B_Mo');
% TET99_B = [];
% TET99_B = xlsread('FossilCoralInput.xls','TET99_B');

% Copy in the name of the coral in these two locations for the one which is
% uncommented out above
% current = TMN;
% samplename = '05TMN';
current = TMI1;
samplename = '11-TM-I1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this line reads in the current timeseries and variable picked
% for the ts calculation
ts = [current(:,1) current(:,variable)];

num = numberpointsperyear - 1;

numbins = round(length(ts)/numberpointsperyear);

% create new interpolated timeseries
% Constructing a cubic spline fit to the data using numknots (# of knots)
% = the number of bins * an arbitrary splinesensitivity input above
numknots = numbins*splinesensitivity;

drilldepth = ts(:,1);
d18Odepth = ts(:,2);
for ww = 1:length(d18Odepth)
    if isnan(d18Odepth(ww)) == 1;
        w(ww) = 0;
    else
        w(ww) = 1;
    end
end
w = w';
sp=spap2(numknots,3,drilldepth,d18Odepth,w);

% gets rid of spurious 0's due to NaN's in original data to make NaN in the
% spline curve for plotting purposes
splined18O = fnval(sp,drilldepth);
for  ww = 1:length(d18Odepth)
    if splined18O(ww) == 0;
        splined18O(ww) = NaN;
    else
        splined18O(ww) = splined18O(ww);
    end
end

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
plot(ts(:,1),ts(:,2),'k-')

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

% %%
% If one manually changes the summers and winters then run from below to
% the end to recalculate the seasonal cycle calculations
%
% Uncomment the lines below if the summers and winters are manually changed
% then run just this cell to recalculate the seasonal cycle calc

% close all
% clear  newdepth scnolag plotnolag sclag plotlag coralinterp coralts time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolating the timeseries to equal monthly spacing using the peaks and
% valleys, or summers and winters picked previously to make 5 datapoints
% between each peak and valley to have 12 points/months in a year

counttime = 1;
% Case 1
if summer(1,1) < winter (1,1) && length(summer) > length(winter)
    for i = 1:(length(summer))
        if i < length(summer)
            newdepth(counttime:counttime+5) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1)-(winter(i,1)-summer(i,1))/6;
            newdepth(counttime+6:counttime+11) = winter(i,1):(summer(i+1,1)-winter(i,1))/6:summer(i+1)-(summer(i+1,1)-winter(i,1))/6;
        elseif i == length(summer)
            xxxx = length(newdepth);
            newdepth(xxxx+1) = summer(length(summer(:,1)),1);
        end
        counttime = counttime + 12;
    end
% Case 2
elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
    for i = 1:(length(summer))
        if i < length(summer)
            newdepth(counttime:counttime+5) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1)-(winter(i,1)-summer(i,1))/6;
            newdepth(counttime+6:counttime+11) = winter(i,1):(summer(i+1,1)-winter(i,1))/6:summer(i+1)-(summer(i+1,1)-winter(i,1))/6;
        elseif i == length(summer)
            newdepth(counttime:counttime+6) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1);
        end
        counttime = counttime + 12;
    end
% Case 3
elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
    for i = 1:(length(winter))
        if i < length(winter)
            newdepth(counttime:counttime+5) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i)-(summer(i,1)-winter(i,1))/6;
            newdepth(counttime+6:counttime+11) = summer(i,1):(winter(i+1,1)-summer(i,1))/6:winter(i+1)-(winter(i+1,1)-summer(i,1))/6;
        elseif i == length(winter)
            xxxx = length(newdepth);
            newdepth(xxxx+1) = winter(length(winter(:,1)),1);
        end
        counttime = counttime + 12;
    end
% Case 4
elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
    for i = 1:(length(winter))
        if i < length(winter)
            newdepth(counttime:counttime+5) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i)-(summer(i,1)-winter(i,1))/6;
            newdepth(counttime+6:counttime+11) = summer(i,1):(winter(i+1,1)-summer(i,1))/6:winter(i+1)-(winter(i+1,1)-summer(i,1))/6;
        elseif i == length(winter)
            newdepth(counttime:counttime+6) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i);
        end
        counttime = counttime + 12;
    end
end

newdepth = newdepth';

% adding in the timestep before the first peak and after the last to
% provide three timesteps for the seasonal averaging
if summer(1,1) < winter(1,1)
    [rb,cb,vb] = find(ts(:,1) == summer(1,1));
    if summer(length(summer),1) > winter(length(winter),1)
        [re,ce,ve] = find(ts(:,1) == summer(length(summer),1));
    elseif summer(length(summer),1) < winter(length(winter),1)
        [re,ce,ve] = find(ts(:,1) == winter(length(winter),1));
    end
elseif summer(1,1) > winter(1,1)
    [rb,cb,vb] = find(ts(:,1) == winter(1,1));
    if summer(length(summer),1) > winter(length(winter),1)
        [re,ce,ve] = find(ts(:,1) == summer(length(summer),1));
    elseif summer(length(summer),1) < winter(length(winter),1)
        [re,ce,ve] = find(ts(:,1) == winter(length(winter),1));
    end
end

coralinterp = interp1(ts(:,1),ts(:,2),newdepth,'spline','extrap');
coralts = [ts(rb-1,2);coralinterp;ts(re+1,2)];
% timestarts on the month before the first seasonal extreme and ends on the
% month after the last seasonal extreme

time = 1:length(coralts);
time = time';

% scnolag = the first summer or winter is difference from the next;
% if the seasons are not the same length, then one is shortened to be the
% same length as the other

ct = 1;
% Case 1
if summer(1,1) < winter (1,1) && length(summer) > length(winter)
    for k = 1:length(summer)-1
        scnolag(k) = (coralts(ct+1,1)) - (coralts(ct+7,1));
        plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
        sclag(k) = (coralts(ct+13,1)) - (coralts(ct+7,1));
        plotlag(k) = (summer(k+1,1) + winter(k,1))/2.;
        ct = ct + 12;
    end
% Case 2
elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
    for k = 1:length(summer)
        scnolag(k) = (coralts(ct+1,1)) - (coralts(ct+7,1));
        plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
        ct = ct + 12;
    end
    ct = 1;
    for k = 1:length(summer)-1
        sclag(k) = (coralts(ct+13,1)) - (coralts(ct+7,1));
        plotlag(k) = (summer(k+1,1) + winter(k,1))/2.;
        ct = ct + 12;
    end
% Case 3    
elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
    for k = 1:length(summer)-1
        scnolag(k) = (coralts(ct+7,1)) - (coralts(ct+1,1));
        plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
        sclag(k) = (coralts(ct+7,1)) - (coralts(ct+13,1));
        plotlag(k) = (summer(k,1) + winter(k+1,1))/2.;
        ct = ct + 12;
    end
% Case 4    
elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
    for k = 1:length(summer)
        scnolag(k) = (coralts(ct+7,1)) - (coralts(ct+1,1));
        plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
        ct = ct + 12;
    end
    ct = 1;
    for k = 1:length(summer)-1
        sclag(k) = (coralts(ct+7,1)) - (coralts(ct+13,1));
        plotlag(k) = (summer(k,1) + winter(k+1,1))/2.;
        ct = ct + 12;
    end
end

scnolag = scnolag';
sclag = sclag';

ssncycnolag = mean(scnolag);
ssncyclag = mean(sclag);

% Figure to check if the first derivative picks summers and winters
% correctly where the original and spline timeseries are also plotted
figure(2)
hold on
set(gca,'Ylim', [8.7 9.2], 'YDir','reverse','Xlim',[0 max(ts(:,1)+10)])
% set(gca,'XLim',[17 45]);
plot(drilldepth,splined18O,'b-')
% plot(ts(:,1),ts(:,2),'k-')
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
set(gca,'YDir','reverse','Xlim',[0 max(time(:,1)+10)])
plot(time,coralts)
hold off

figure(4)
subplot(3,1,1)
hold on
plot(ts(:,1),ts(:,2),'k-')
plot(ts(:,1),splined18O,'b-')
plot(winter(:,1),winter(:,2),'bo')
plot(summer(:,1),summer(:,2),'ro')
set(gca,'Xlim',[0 max(ts(:,1)+10)],'YDir','reverse','Ylim', [8.7 9.2],'FontName','Myriad Pro','fontsize',16);
% set(gca,'XLim',[17 45],'YDir','reverse','FontName','Myriad Pro','fontsize',16);
%xlabel('Sample Number')
xlabel('Depth (cm)')
%ylabel({'Coral \delta^{18}O';'(permil)'})
ylabel({'Coral Sr/Ca';'(mmol/mol)'})
hold off
subplot(3,1,2)
hold on
bar(plotnolag,-1*scnolag,1,'k')
line([0 max(ts(:,1)+10)],[-1*ssncycnolag -1*ssncycnolag])
set(gca,'Xlim',[0 max(ts(:,1)+10)])
% set(gca,'XLim',[17 45]);
set(gca,'FontName','Myriad Pro','fontsize',16);
ylabel({'Seasonal Cycle';'Amplitude'})
xlabel('Sample Number')
hold off
subplot(3,1,3)
hold on
bar(plotlag,-1*sclag,1,'g')
line([0 max(ts(:,1)+10)],[-1*ssncyclag -1*ssncyclag])
set(gca,'Xlim',[0 max(ts(:,1)+10)])
% set(gca,'XLim',[17 45]);
set(gca,'FontName','Myriad Pro','fontsize',16);
ylabel({'Seasonal Cycle';'Amplitude'})
xlabel('Sample Number')
hold off
% print('-cmyk','-dtiff','-r300','SCpickerexample')

% Doing a bootstrap analysis to determine the 95% confidence intervals for
% the mean of the computed climatologies
scmeannolag = @(scnolag)mean(scnolag);
scmeanlag = @(sclag)mean(sclag);

nolagbootci = bootci(1000,scmeannolag,scnolag);
lagbootci = bootci(1000,scmeanlag,sclag);

avgsc = -1*(ssncycnolag + ssncyclag)/2.;
avgerror = -1*((nolagbootci(1,1) - ssncycnolag) + (ssncycnolag - nolagbootci(2,1)) + (lagbootci(1,1) - ssncyclag) + (ssncyclag - lagbootci(2,1)))/4.;
disp(['The number of complete years in the record ',samplename,' is ',num2str(numberofyearsinrecord,3),' years'])
if variable == 2;
    temp = avgsc/0.22;
    temperror = avgerror/0.22;
    disp(['The d18O seasonal cycle = ',num2str(avgsc,2),' +/- ',num2str(avgerror,2)])
    disp(['In degrees, the seasonal cycle = ',num2str(temp,2),' +/- ',num2str(temperror,2)])
elseif variable == 3;
    temp = avgsc/0.058;
    temperror = avgerror/0.058;
    disp(['The Sr/Ca seasonal cycle = ',num2str(avgsc,2),' +/- ',num2str(avgerror,2)])
    disp(['In degrees, the seasonal cycle = ',num2str(temp,2),' +/- ',num2str(temperror,2)])
elseif variable == 4;
    disp(['The d13C seasonal cycle = ',num2str(avgsc,2),' +/- ',num2str(avgerror,2)])
end

%======================================
% Coral Timeseries Input for Event Picker
%======================================

coralts=[time/12.+(11/12) coralts(:,1)]; 

% %======================================
% % Coral geochemical timeseries and climatology
% %======================================
% %Output the coral Sr/Ca timeseries 
% %variable = 3 (Sr/Ca)
% %variable = 2 (d18O)
% 
% if variable == 3; 
%     SrCatimeseries = coralts;
%     clima = climatology(variable,SrCatimeseries(:,2));
% elseif variable == 2;
%     d18Otimeseries = coralts;
%     clima = climatology(variable,d18Otimeseries(:,2));
% elseif variable == 3;
%     d13Ctimeseries = coralts;
%     clima = climatology(variable,d13Ctimeseries(:,2));
% end
% 
% %======================================
% % Coral climatology
% %======================================
% clima = climatology(variable,SrCatimeseries(:,2));