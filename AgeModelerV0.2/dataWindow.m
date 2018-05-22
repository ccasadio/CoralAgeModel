classdef dataWindow < handle
    
    properties
        container;
        value;
        text;
    end
    
    events
        dataUpdated;
    end
    
    methods
        function obj = dataWindow(parent, position, name)
            
            obj.container = uipanel(parent,'Title', name, 'Position',position);
            
            obj.text = uicontrol(obj.container,'Style','text',...
                'FontUnits','normalized',...
                'FontSize',.4,...
                'Units','normalized',...
                'Position',[.1 .1 .8 .8]);
            
        end
        
        function setValue(obj,newValue)
            obj.value = newValue;
            obj.text.String = num2str(newValue);
            notify(obj,'dataUpdated');
        end
        
        function l = addDataListener(obj,funHandle)
            l = addlistener(obj,'dataUpdated',funHandle);
        end
    end
end

