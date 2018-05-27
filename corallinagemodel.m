function [out] = corallinagemodel( inputData )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[pow,w] = pwelch(inputData(:,2));
[~, loc] = findpeaks(pow, w, 'SortStr', 'descend');%'NPeaks', 1,
loc = 2*pi./loc;
loc = loc((loc < 36) & (loc > 3));
ppy = loc(1);

out = linearAgeModel(inputData,ppy,2.7);
end

function [data] = linearAgeModel(inputData, numberpointsperyear, splinesensitivity)

    numbins = round(length(inputData)/numberpointsperyear)

    % create new interpolated timeseries
    % Constructing a cubic spline fit to the data using numknots (# of knots)
    % = the number of bins * an arbitrary splinesensitivity input above
    numknots = numbins*splinesensitivity

    for i = 1:length(inputData(:,2))
        if isnan(inputData(i,2)) == 1;
            w(i) = 0;
        else
            w(i) = 1;
        end
    end
    w = w';
    sp=spap2(numknots,3,inputData(:,1),inputData(:,2),w);

    % gets rid of spurious 0's due to NaN's in original data to make NaN in the
    % spline curve for plotting purposes
    splined18O = fnval(sp,inputData(:,1));
    for  i = 1:length(inputData(:,2))
        if splined18O(i) == 0;
            splined18O(i) = NaN;
        else
            splined18O(i) = splined18O(i);
        end
    end

    % Calculating the first derivative of the spline fit
    firstder = fnder(sp,1);
    fderval = fnval(firstder,inputData(:,1));

    % Logging the inflection points where the first derivative = 0 (sign
    % change)
    counter = 1;
    for i = 1:length(fderval)-1
        if fderval(i,1) > 0 && fderval(i+1,1) < 0
            if isnan(inputData(i,2)) == 1
                inflections(counter,:) = [inputData(i,1) splined18O(i)];
            else
                inflections(counter,:) = inputData(i,:);
            end
            counter = counter + 1;
        elseif fderval(i,1) < 0 && fderval(i+1,1) > 0
            if isnan(inputData(i,2)) == 1
                inflections(counter,:) = [inputData(i,1) splined18O(i)];
            else
                inflections(counter,:) = inputData(i,:);
            end
            counter = counter + 1;
        end
    end

    % computing mean "distance" for a year. This number should be close to the
    % number input by the user in the beginning of the program
    for i = 1:length(inflections)-1
        meanyears(i,1) = inflections(i+1,1)-inflections(i,1);
    end
    meandistance = mean(meanyears);

    % Number of years in the record which is half of the number of inflections
    % points
    numberofyearsinrecord = length(inflections)/2;

    % computing seasonal cycle
    % compiling all of the winters, which is the max value in the
    % surrounding five data points from the inflection point
    % If data gaps are filled using NaN, then the splined18O value is used 
    counter = 1;
    for i = 2:length(fderval)-1
        if fderval(i,1) > 0 && fderval(i+1,1) < 0
            if isnan(inputData(i-2,2))||isnan(inputData(i-1,2))||isnan(inputData(i,2))||isnan(inputData(i+1,2))||isnan(inputData(i+2,2)) == 1;
                if isnan(inputData(i-2,2)) == 1
                    newminus2 = [inputData(i-2,1) splined18O(i-2,:)];
                else
                    newminus2 = inputData(i-2,:);
                end
                if isnan(inputData(i-1,2)) == 1
                    newminus1 = [inputData(i-1,1) splined18O(i-1,:)];
                else
                    newminus1 = inputData(i-1,:);
                end
                if isnan(inputData(i,2)) == 1
                    newone = [inputData(i,1) splined18O(i,:)];
                else
                    newone = inputData(i,:);
                end
                if isnan(inputData(i+1,2)) == 1
                    newplus1 = [inputData(i+1,1) splined18O(i+1,:)];
                else
                    newplus1 = inputData(i+1,:);
                end
                if isnan(inputData(i+1,2)) == 1
                    newplus2 = [inputData(i+2,1) splined18O(i+2,:)];
                else
                    newplus2 = inputData(i+2,:);
                end
                winterblock = [newminus2;newminus1;newone;newplus1;newplus2];
                [r,~,~] = find(winterblock == max(winterblock(:,2)));
                winter(counter,:) = winterblock(r(1,1),:);
            else
                winterblock = [inputData(i-2,:);inputData(i-1,:);inputData(i,:);inputData(i+1,:);inputData(i+2,:)];
                [r,~,~] = find(winterblock == max(winterblock(:,2)));
                winter(counter,:) = winterblock(r(1,1),:);
            end
            counter = counter + 1;
        end
    end

    % compiling all of summers, which is the min value in the
    % surrounding five data points from the inflection point
    % If data gaps are filled using NaN, then the splined18O value is used
    counter = 1;
    for i = 2:length(fderval)-1
        if fderval(i,1) < 0 && fderval(i+1,1) > 0
            if isnan(inputData(i-2,2))||isnan(inputData(i-1,2))||isnan(inputData(i,2))||isnan(inputData(i+1,2))||isnan(inputData(i+2,2)) == 1;
                if isnan(inputData(i-1,2)) == 1
                    newminus2 = [inputData(i-2,1) splined18O(i-2,:)];
                else
                    newminus2 = inputData(i-2,:);
                end
                if isnan(inputData(i-1,2)) == 1
                    newminus1 = [inputData(i-1,1) splined18O(i-1,:)];
                else
                    newminus1 = inputData(i-1,:);
                end
                if isnan(inputData(i,2)) == 1
                    newone = [inputData(i,1) splined18O(i,:)];
                else
                    newone = inputData(i,:);
                end
                if isnan(inputData(i+1,2)) == 1
                    newplus1 = [inputData(i+1,1) splined18O(i+1,:)];
                else
                    newplus1 = inputData(i+1,:);
                end
                if isnan(inputData(i+2,2)) == 1
                    newplus2 = [inputData(i+2,1) splined18O(i+2,:)];
                else
                    newplus2 = inputData(i+2,:);
                end
                summerblock = [newminus2;newminus1;newone;newplus1;newplus2];
                [r,~,~] = find(summerblock == min(summerblock(:,2)));
                summer(counter,:) = summerblock(r(1,1),:);
            else
                summerblock = [inputData(i-2,:);inputData(i-1,:);inputData(i,:);inputData(i+1,:);inputData(i+2,:)];
                [r,~,~] = find(summerblock == min(summerblock(:,2)));
                summer(counter,:) = summerblock(r(1,1),:);
            end
            counter = counter + 1;
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

    count = 1;
    % Case 1
    if summer(1,1) < winter (1,1) && length(summer) > length(winter)
        for i = 1:(length(summer))
            if i < length(summer)
                newdepth(count:count+5) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1)-(winter(i,1)-summer(i,1))/6;
                newdepth(count+6:count+11) = winter(i,1):(summer(i+1,1)-winter(i,1))/6:summer(i+1)-(summer(i+1,1)-winter(i,1))/6;
            elseif i == length(summer)
                xxxx = length(newdepth);
                newdepth(xxxx+1) = summer(length(summer(:,1)),1);
            end
            count = count + 12;
        end
    % Case 2
    elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
        for i = 1:(length(summer))
            if i < length(summer)
                newdepth(count:count+5) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1)-(winter(i,1)-summer(i,1))/6;
                newdepth(count+6:count+11) = winter(i,1):(summer(i+1,1)-winter(i,1))/6:summer(i+1)-(summer(i+1,1)-winter(i,1))/6;
            elseif i == length(summer)
                newdepth(count:count+6) = summer(i,1):(winter(i,1)-summer(i,1))/6:winter(i,1);
            end
            count = count + 12;
        end
    % Case 3
    elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
        for i = 1:(length(winter))
            if i < length(winter)
                newdepth(count:count+5) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i)-(summer(i,1)-winter(i,1))/6;
                newdepth(count+6:count+11) = summer(i,1):(winter(i+1,1)-summer(i,1))/6:winter(i+1)-(winter(i+1,1)-summer(i,1))/6;
            elseif i == length(winter)
                xxxx = length(newdepth);
                newdepth(xxxx+1) = winter(length(winter(:,1)),1);
            end
            count = count + 12;
        end
    % Case 4
    elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
        for i = 1:(length(winter))
            if i < length(winter)
                newdepth(count:count+5) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i)-(summer(i,1)-winter(i,1))/6;
                newdepth(count+6:count+11) = summer(i,1):(winter(i+1,1)-summer(i,1))/6:winter(i+1)-(winter(i+1,1)-summer(i,1))/6;
            elseif i == length(winter)
                newdepth(count:count+6) = winter(i,1):(summer(i,1)-winter(i,1))/6:summer(i);
            end
            count = count + 12;
        end
    end

    newdepth = newdepth';

    % adding in the timestep before the first peak and after the last to
    % provide three timesteps for the seasonal averaging
    if summer(1,1) < winter(1,1)
        [rb,~,~] = find(inputData(:,1) == summer(1,1));
        if summer(length(summer),1) > winter(length(winter),1)
            [re,~,~] = find(inputData(:,1) == summer(length(summer),1));
        elseif summer(length(summer),1) < winter(length(winter),1)
            [re,~,~] = find(inputData(:,1) == winter(length(winter),1));
        end
    elseif summer(1,1) > winter(1,1)
        [rb,~,~] = find(inputData(:,1) == winter(1,1));
        if summer(length(summer),1) > winter(length(winter),1)
            [re,~,~] = find(inputData(:,1) == summer(length(summer),1));
        elseif summer(length(summer),1) < winter(length(winter),1)
            [re,~,~] = find(inputData(:,1) == winter(length(winter),1));
        end
    end

    coralinterp = interp1(inputData(:,1),inputData(:,2),newdepth,'spline','extrap');
    data = [inputData(rb-1,2);coralinterp;inputData(re+1,2)];
    % timestarts on the month before the first seasonal extreme and ends on the
    % month after the last seasonal extreme

    time = 1:length(data);
    time = time';

    % scnolag = the first summer or winter is difference from the next;
    % if the seasons are not the same length, then one is shortened to be the
    % same length as the other

    ct = 1;
    % Case 1
    if summer(1,1) < winter (1,1) && length(summer) > length(winter)
        for k = 1:length(summer)-1
            scnolag(k) = (data(ct+1,1)) - (data(ct+7,1));
            plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
            sclag(k) = (data(ct+13,1)) - (data(ct+7,1));
            plotlag(k) = (summer(k+1,1) + winter(k,1))/2.;
            ct = ct + 12;
        end
    % Case 2
    elseif summer(1,1) < winter (1,1) && length(summer) == length(winter)
        for k = 1:length(summer)
            scnolag(k) = (data(ct+1,1)) - (data(ct+7,1));
            plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
            ct = ct + 12;
        end
        ct = 1;
        for k = 1:length(summer)-1
            sclag(k) = (data(ct+13,1)) - (data(ct+7,1));
            plotlag(k) = (summer(k+1,1) + winter(k,1))/2.;
            ct = ct + 12;
        end
    % Case 3    
    elseif summer(1,1) > winter (1,1) && length(summer) < length(winter)
        for k = 1:length(summer)-1
            scnolag(k) = (data(ct+7,1)) - (data(ct+1,1));
            plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
            sclag(k) = (data(ct+7,1)) - (data(ct+13,1));
            plotlag(k) = (summer(k,1) + winter(k+1,1))/2.;
            ct = ct + 12;
        end
    % Case 4    
    elseif summer(1,1) > winter (1,1) && length(summer) == length(winter)
        for k = 1:length(summer)
            scnolag(k) = (data(ct+7,1)) - (data(ct+1,1));
            plotnolag(k) = (summer(k,1) + winter(k,1))/2.;
            ct = ct + 12;
        end
        ct = 1;
        for k = 1:length(summer)-1
            sclag(k) = (data(ct+7,1)) - (data(ct+13,1));
            plotlag(k) = (summer(k,1) + winter(k+1,1))/2.;
            ct = ct + 12;
        end
    end

    scnolag = scnolag';
    sclag = sclag';

    ssncycnolag = mean(scnolag);
    ssncyclag = mean(sclag);



    %======================================
    % Coral Timeseries Input for Event Picker
    %======================================

    data=[time/12 data(:,1)];
end