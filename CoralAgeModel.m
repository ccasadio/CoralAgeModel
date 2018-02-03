function varargout = CoralAgeModel(varargin)
% CORALAGEMODEL MATLAB code for CoralAgeModel.fig
%      CORALAGEMODEL, by itself, creates a new CORALAGEMODEL or raises the existing
%      singleton*.
%
%      H = CORALAGEMODEL returns the handle to a new CORALAGEMODEL or the handle to
%      the existing singleton*.
%
%      CORALAGEMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORALAGEMODEL.M with the given input arguments.
%
%      CORALAGEMODEL('Property','Value',...) creates a new CORALAGEMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CoralAgeModel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CoralAgeModel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CoralAgeModel

% Last Modified by GUIDE v2.5 17-Jan-2018 12:28:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoralAgeModel_OpeningFcn, ...
                   'gui_OutputFcn',  @CoralAgeModel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

assert(nargin > 0, 'CoralAgeModel requires input data.');

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before CoralAgeModel is made visible.
function CoralAgeModel_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CoralAgeModel (see VARARGIN)

% Choose default command line output for CoralAgeModel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.spline_sensitivity_text.String = num2str(3);
handles.points_per_year_text.String = num2str(12);
handles.data = varargin{1};

% gets rid of spurious 0's due to NaN's in original data to make NaN in the
% spline curve for plotting purposes
handles.data(handles.data(:,2) == 0, 2) = NaN;

handles.tgroup = uitabgroup('Parent', handles.panel);
handles.tab1 = axes('Parent',uitab('Parent', handles.tgroup, 'Title', 'Spline Fit'));
handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Inflection Points');
handles.tab3 = uitab('Parent', handles.tgroup, 'Title', 'Seasonal Cycle');
handles.tab4 = uitab('Parent', handles.tgroup, 'Title', 'Peak to Peak');
handles.tab5 = uitab('Parent', handles.tgroup, 'Title', 'Descriptive Stats');

guidata(hObject, handles); 

splineFit(handles);


% --- Outputs from this function are returned to the command line.
function varargout = CoralAgeModel_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in increase_spline_sensitivity_button.
function increase_spline_sensitivity_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to increase_spline_sensitivity_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
splineSensitivity = str2double(handles.spline_sensitivity_text.String) + .1;
handles.spline_sensitivity_text.String = num2str(splineSensitivity);
splineFit(handles);

% --- Executes on button press in decrease_spline_sensitivity_button.
function decrease_spline_sensitivity_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to decrease_spline_sensitivity_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
splineSensitivity = str2double(handles.spline_sensitivity_text.String) - .1;
handles.spline_sensitivity_text.String = num2str(splineSensitivity);
splineFit(handles);

% --- Executes on button press in increase_points_per_year_button.
function increase_points_per_year_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to increase_points_per_year_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pointsPerYear = str2double(handles.points_per_year_text.String) + 1;
handles.points_per_year_text.String = num2str(pointsPerYear);
splineFit(handles);

% --- Executes on button press in decrease_points_per_year_button.
function decrease_points_per_year_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to decrease_points_per_year_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pointsPerYear = str2double(handles.points_per_year_text.String) - 1;
handles.points_per_year_text.String = num2str(pointsPerYear);
splineFit(handles);

function splineFit(handles)

data = handles.data;
splineSensitivity = str2double(handles.spline_sensitivity_text.String);
pointsPerYear = str2double(handles.points_per_year_text.String);

%% Create new interpolated timeseries
%  Constructing a cubic spline fit to the data using numknots (# of knots)
%  = the number of bins * an arbitrary splinesensitivity input above
numknots = (length(data)/pointsPerYear)*splineSensitivity;

warning('off','SPLINES:CHCKXYWP:NaNs');
sp = spap2(numknots,3,data(:,1),data(:,2));
splined18O = fnval(sp,data(:,1));

% Calculating the first derivative of the spline fit
firstder = fnder(sp,1);
fderval = fnval(firstder,data(:,1));

%inflection points of the spline fit to the original data
inflectionPointsSpline = fnzeros(firstder)';
inflectionPointsSpline(:,2) = fnval(sp, inflectionPointsSpline(:,1));

%method 1 to find data inflection points:
%   use x value of spline inflection points to find the closest real data x
%   value. Doesn't garuntee the largest local y value.
[ ~, inflectionPointIndicies] = arrayfun(@(x)min(abs(x - data(:,1))), inflectionPointsSpline(:,1));
inflectionPointsData1 = [data(inflectionPointIndicies,1), data(inflectionPointIndicies,2)];

%method 2 to find data inflection points:
%   use findpeaks to search y values of original data to find minima and
%   maxima.
meanData = nanmean(data(:,2));
[ ~, maxInflectionPointIndicies] = findpeaks(data(:,2),'MinPeakDistance', pointsPerYear*3/4, 'MinPeakHeight', meanData);
[ ~, minInflectionPointIndicies] = findpeaks(data(:,2) * -1,'MinPeakDistance', pointsPerYear*3/4, 'MinPeakHeight', meanData * -1);
inflectionPointIndicies = sort([maxInflectionPointIndicies; minInflectionPointIndicies]);
inflectionPointsData2 = data(inflectionPointIndicies,:);

%% Calculate Metrics
meanYearsSpline = 2 * mean(diff(inflectionPointsSpline(:,1)));
meanYearsData = 2 * mean(diff(inflectionPointsData1(:,1)));
meanYearsData2 = 2 * mean(diff(inflectionPointsData2(:,1)));


% Logging the inflection points where the first derivative = 0 (sign
% change)
judsSummerPick = fderval(1:(end-1),1) > 0 & fderval(2:end,1) < 0;
judsWinterPick = fderval(1:(end-1),1) < 0 & fderval(2:end,1) > 0;
inflections = data(judsSummerPick | judsWinterPick, :);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                   Updating UI                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Metrics
handles.year_length_mean_spline.String = num2str(meanYearsSpline);
handles.year_length_mean_peaks.String = num2str(meanYearsData2);

handles.number_of_calculated_years_peaks.String = num2str((range(data(:,1)))/meanYearsData2);
handles.number_of_calculated_years_spline.String = num2str(range(data(:,1))/meanYearsSpline);

seasonalCyclePeaks = mean(abs(diff(inflectionPointsData2(:,2))));
seasonalCycleSpline = mean(abs(diff(inflectionPointsSpline(:,2))));

handles.seasonal_cycle_peaks.String = num2str(seasonalCyclePeaks);
handles.seasonal_cycle_spline.String = num2str(seasonalCycleSpline);

%%
% Spline fit to data
plot(handles.tab1,data(:,1),data(:,2),'r-o','linewidth',.2,'markersize',3);
hold(handles.tab1, 'on');
plot(handles.tab1,data(:,1),splined18O,'k-','linewidth',1.4);
title(handles.tab1, 'Spline Fit to Data');
hold(handles.tab1, 'off');

% Peaks on data
ax(1) = subplot(3,1,1,'parent',handles.tab2);
plot(ax(1),data(:,1),data(:,2),'b-');
hold(ax(1), 'on');
plot(ax(1),inflectionPointsSpline(:,1),inflectionPointsSpline(:,2),'g*');
title(ax(1), 'Spline Inflection Points');
hold(ax(1), 'off');
ax(2) = subplot(3,1,2,'parent',handles.tab2);
plot(ax(2),data(:,1),data(:,2),'b-');
hold(ax(2), 'on');
plot(ax(2),inflectionPointsData2(:,1),inflectionPointsData2(:,2),'g*');
title(ax(2), '"Peaks" Inflection Points');
hold(ax(2), 'off');
ax(3) = subplot(3,1,3,'parent',handles.tab2);
plot(ax(3),data(:,1),data(:,2),'b-');
hold(ax(3), 'on');
plot(ax(3),inflectionPointsData1(:,1),inflectionPointsData1(:,2),'g*');
title(ax(3), 'Closest Data to Spline Inflection');
hold(ax(3), 'off');

linkaxes(ax);
clear ax;

% Seasonal Cycle
ax(1) = subplot(4,1,1,'parent',handles.tab3);
plot(ax(1),data(:,1),data(:,2),'b-');
title(ax(1), 'Data');
ax(3) = subplot(4,1,2,'parent',handles.tab3);
bar(ax(3),(inflectionPointsData2(2:end,1)+inflectionPointsData2(1:(end-1),1))/2,diff(inflectionPointsData2(:,2)));
title(ax(3), '"Peaks" Seasonal Cycle');
try
    ax(2) = subplot(4,1,3,'parent',handles.tab3);
    plot(ax(2),data(:,1),splined18O,'b-');
    title(ax(2), 'Spline');
    ax(4) = subplot(4,1,4,'parent',handles.tab3);
    bar(ax(4),(inflectionPointsSpline(2:end,1)+inflectionPointsSpline(1:(end-1),1))/2,diff(inflectionPointsSpline(:,2)));
    title(ax(4), 'Spline Seasonal Cycle');
catch e
    disp(e.message);
    [~, i] = unique(inflectionPointsSpline(2:end,1));
    disp(inflectionPointsSpline(not(ismember(1:numel(inflectionPointsSpline(2:end,1)),i)),1));
end

linkaxes(ax(1:2));
linkaxes(ax(3:4));
clear ax;


% Peak to Peak
ax(1) = subplot(3,1,1,'parent',handles.tab4);
plot(ax(1),data(:,1),data(:,2),'b-');
title(ax(1), 'Data');
ax(2) = subplot(3,1,2,'parent',handles.tab4);
plot(ax(2),inflectionPointsData1(:,1),inflectionPointsData1(:,2),'b-');
title(ax(2), 'Spline Peak to Peak');
ax(3) = subplot(3,1,3,'parent',handles.tab4);
plot(ax(3),inflectionPointsData2(:,1),inflectionPointsData2(:,2),'b-');
title(ax(3), '"Peaks" Peak to Peak');

linkaxes(ax);
clear ax;

% Desctiptive stats
ax(1) = subplot(3,2,1,'parent',handles.tab5);
ax(2) = subplot(3,2,2,'parent',handles.tab5);
ax(3) = subplot(3,2,3,'parent',handles.tab5);
ax(4) = subplot(3,2,4,'parent',handles.tab5);
histogram(ax(1),data(:,2), 'Normalization', 'probability');
histogram(ax(2),splined18O, 'Normalization', 'probability');
histogram(ax(3),abs(diff(inflectionPointsData2(:,2))), 'Normalization', 'probability');
histogram(ax(4),abs(diff(inflectionPointsSpline(:,2))), 'Normalization', 'probability');

ax(1).YLim = [0 .4];
ax(2).YLim = [0 .4];
ax(3).YLim = [0 .4];
ax(4).YLim = [0 .4];

clear ax;