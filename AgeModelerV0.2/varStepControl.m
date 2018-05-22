classdef varStepControl < handle
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        value;
        delta;
        container;
    end
    
    events
        dataUpdated;
    end
    
    methods
        function obj = varStepControl(parent, position, name, startVal, startDelta, deltaDelta)
            obj.value = startVal;
            obj.delta = startDelta;
            obj.container = uipanel(parent,'Title', name, 'Position',position);
            c1 = stepControl(obj.container,[.05 0 .9 .5],'step',startVal,startDelta);
            c2 = stepControl(obj.container,[.05 .5 .9 .5],[char(916),' step'],startDelta,deltaDelta);
            c2.addDataListener(@(o,~) c1.setDelta(o.value));
            
            function link(o,~)
                obj.value = o.value;
                notify(obj,'dataUpdated');
            end
            
            c1.addDataListener(@link);
        end
        
        function l = addDataListener(obj,funHandle)
            l = addlistener(obj,'dataUpdated',funHandle);
        end
    end
    
end

