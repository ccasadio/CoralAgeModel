classdef spAgeModelFitter < handle
    
    properties
        fig;
        timeSeries;
        ageModel;
        criticalPoints;
        tsPowerSpectrum;
        amPowerSpectrum;
        splineSensitivity = 2.7;
        pointsPerYear = 12;
    end
    
    events
        tsUpdated;
        amUpdated;
        ssUpdated;
        ppyUpdated;
    end
    
    methods
        
        function obj = spAgeModelFitter(ts)
            %create main figure
            obj.fig = figure;
            
            controlHeight = .23;
            
            tg = uitabgroup(obj.fig, 'Position', [0 controlHeight 1 (1-controlHeight)]);
            
            %% Tab 1 - Basic Age Model View 
            tabs(1) = uitab(tg, 'Title', 'Age Model');
            
            axesPosition2x1 = {[.05 0.575 .9 .375], [.05 0.1 .9 .375]};
            
            %create axes for plots of timeSeries and ageModel
            axTS = axes('Parent',tabs(1),'Units','normalized','Position',axesPosition2x1{1});
            axAM = axes('Parent',tabs(1),'Units','normalized','Position',axesPosition2x1{2});
            
            function plotTS(~,~)
                plot(axTS,obj.timeSeries(:,1),obj.timeSeries(:,2));
                xlim(axTS,[min(obj.timeSeries(:,1)) max(obj.timeSeries(:,1))]);
                title(axTS,'Data');
            end
            
            %when timeSeries is changed, replot it
            obj.addTSListener(@plotTS);

            function plotAM(~,~)
                plot(axAM,obj.ageModel(:,1),obj.ageModel(:,2));
                xlim(axAM,[min(obj.ageModel(:,1)) max(obj.ageModel(:,1))]);
                title(axAM,'Age Model');
            end
            
            %when ageModel is changed, replot it
            obj.addAMListener(@plotAM);

            %recalculate ageModel using values contained in object
            function updateAM(~,~)
                [am, cp] = spAgeModel(obj.timeSeries, obj.pointsPerYear, obj.splineSensitivity);
                obj.criticalPoints = cp;
                obj.setAM(am);
            end
            
            %recalculate ageModel any time a needed variable is updated
            obj.addTSListener(@updateAM);
            obj.addSSListener(@updateAM);
            obj.addPPYListener(@updateAM);
            
            %% Tab 2 - Frequency Analysis 
            tabs(2) = uitab(tg, 'Title', 'Frequency Comp');
            
            %create axes for plots of timeSeries and ageModel power spcetrum
            axTSPowSp = axes('Parent',tabs(2),'Units','normalized','Position',axesPosition2x1{1});
            axAMPowSp = axes('Parent',tabs(2),'Units','normalized','Position',axesPosition2x1{2});
            
            function adjustAxes()
                xlim(axTSPowSp,[min([obj.tsPowerSpectrum(:,1); obj.amPowerSpectrum(:,1)]), max([obj.tsPowerSpectrum(:,1); obj.amPowerSpectrum(:,1)])]);
                ylim(axTSPowSp,[min([obj.tsPowerSpectrum(:,2); obj.amPowerSpectrum(:,2)]), max([obj.tsPowerSpectrum(:,2); obj.amPowerSpectrum(:,2)])]);
                
                xlim(axAMPowSp,[min([obj.tsPowerSpectrum(:,1); obj.amPowerSpectrum(:,1)]), max([obj.tsPowerSpectrum(:,1); obj.amPowerSpectrum(:,1)])]);
                ylim(axAMPowSp,[min([obj.tsPowerSpectrum(:,2); obj.amPowerSpectrum(:,2)]), max([obj.tsPowerSpectrum(:,2); obj.amPowerSpectrum(:,2)])]);
            end
            
            function plotTSPowSp(~,~)
                [pxx,f] = pmtm(obj.timeSeries(:,2));
                obj.tsPowerSpectrum = [2*pi./(obj.pointsPerYear*f), pxx];
                plot(axTSPowSp,obj.tsPowerSpectrum(:,1),obj.tsPowerSpectrum(:,2));
                set(axTSPowSp, 'YScale', 'log', 'XScale', 'log', 'XDir', 'reverse');
                title(axTSPowSp,'Data');
                xlabel(axTSPowSp,'Period');
                ylabel(axTSPowSp,'Power Spectral Density');
                %adjustAxes();
            end
            obj.addTSListener(@plotTSPowSp);
            
            function plotAMPowSp(~,~)    
                [pxx,f] = pmtm(obj.ageModel(:,2));
                obj.amPowerSpectrum = [2*pi./(obj.pointsPerYear*f), pxx];
                plot(axAMPowSp,obj.amPowerSpectrum(:,1),obj.amPowerSpectrum(:,2));
                set(axAMPowSp, 'YScale', 'log', 'XScale', 'log', 'XDir', 'reverse');
                title(axAMPowSp,'Age Model');
                xlabel(axAMPowSp,'Period');
                ylabel(axAMPowSp,'Power Spectral Density');
                adjustAxes();
            end
            obj.addAMListener(@plotAMPowSp);
            
            
             %% Tab 3 - Autocorrelation Analysis 
            tabs(3) = uitab(tg, 'Title', 'Autocorr Comp');
            
            %create axes for plots of timeSeries and ageModel autocorr
            axTSACorr = axes('Parent',tabs(3),'Units','normalized','Position',axesPosition2x1{1});
            axAMACorr = axes('Parent',tabs(3),'Units','normalized','Position',axesPosition2x1{2});
            
            function plotAC(ax,data)
                [acf,lags,bounds] = autocorr(data);
                lineHandles = stem(ax,lags,acf,'filled','r-o');
                hold(ax,'on');
                set(lineHandles(1),'MarkerSize',4)
                grid(ax,'on')
                xlabel(ax,'Lag')
                ylabel(ax,'Sample Autocorrelation')

                plot(ax,[0 0; 20 20],[bounds([1 1]) bounds([2 2])],'-b');
                plot(ax,[0 20],[0 0],'-k');
                
                hold(ax,'off');
            end
            
            function plotTSAC(~,~)
                %  Plot the sample ACF:
                plotAC(axTSACorr,obj.timeSeries(:,2));
                title(axTSACorr,'Data');
            end
            
            %when timeSeries is changed, replot autocorr of it
            obj.addTSListener(@plotTSAC);  
            
            function plotAMAC(~,~)
                %  Plot the sample ACF:
                plotAC(axAMACorr,obj.ageModel(:,2));
                title(axAMACorr,'Age Model');
            end
            
            %when ageModel is changed, replot autocorr of it
            obj.addAMListener(@plotAMAC);
            
            %% Control Section
            controlSpace = uipanel(obj.fig, 'Position', [0 0 1 controlHeight]);
            
            %create a variable step control for splineSensitivity
            vs = varStepControl(controlSpace,[0 0 .2 1],'Spline Sensitivity',obj.splineSensitivity,.1,.01);
            vs.addDataListener(@(o,~) obj.setSS(o.value));
            
            %create a step control for pointsPerYear
            vs = stepControl(controlSpace,[.225 0 .2 .5],'Points Per Year',obj.pointsPerYear,1);
            vs.addDataListener(@(o,~) obj.setPPY(o.value));
            
            %create a data window for the data mean
            dw = dataWindow(controlSpace,[.45 0 .2 .5],'mean ageModel');
            function calc(~,~)
                dw.setValue(mean(obj.ageModel(:,2)));
            end
            obj.addAMListener(@calc);
            
            
            obj.setTS(ts);
        end
        
        function setTS(obj,ts)
            obj.timeSeries = ts;
            notify(obj,'tsUpdated');
        end

        function setAM(obj,am)
            obj.ageModel = am;
            notify(obj,'amUpdated');
        end

        function setSS(obj,ss)
            obj.splineSensitivity = ss;
            notify(obj,'ssUpdated');
        end

        function setPPY(obj,ppy)
            obj.pointsPerYear = ppy;
            notify(obj,'ppyUpdated');
        end
        
        function l = addTSListener(obj,funHandle)
            l = addlistener(obj,'tsUpdated',funHandle);
        end        
        
        function l = addAMListener(obj,funHandle)
            l = addlistener(obj,'amUpdated',funHandle);
        end     
        
        function l = addSSListener(obj,funHandle)
            l = addlistener(obj,'ssUpdated',funHandle);
        end     
        
        function l = addPPYListener(obj,funHandle)
            l = addlistener(obj,'ppyUpdated',funHandle);
        end
    end
    
end

