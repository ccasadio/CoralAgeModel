classdef stepControl < handle
    
    properties
        container;
        value;
        delta;
    end
    
    events
        dataUpdated;
    end
    
    methods
        function obj = stepControl(parent, position, name, startVal, delta)
            obj.value = startVal;
            obj.delta = delta;
            obj.container = uipanel(parent,'Title', name, 'Position',position);
            
            function inc(~, ~, ~)
                obj.value = obj.value + obj.delta;
                obj.container.Children(2).String = num2str(obj.value);
                notify(obj, 'dataUpdated');
            end
            
            function dec(~, ~, ~)
                obj.value = obj.value - obj.delta;
                obj.container.Children(2).String = num2str(obj.value);
                notify(obj, 'dataUpdated');
            end
            
            uicontrol(obj.container,'Style','pushbutton','String',char(9664),...
                'FontUnits','normalized',...
                'Units','normalized',...
                'Position',[.075 .15 .25 .7],...
                'Callback',@dec);
            
            uicontrol(obj.container,'Style','text','String',num2str(startVal),...
                'FontUnits','normalized',...
                'FontSize',.4,...
                'Units','normalized',...
                'Position',[.375 .15 .25 .7]);
            
            uicontrol(obj.container,'Style','pushbutton','String',char(9654),...
                'FontUnits','normalized',...
                'Units','normalized',...
                'Position',[.675 .15 .25 .7],...
                'Callback',@inc);
            
        end
        
        function setDelta(obj,newDel)
            obj.delta = newDel;
        end
        
        function l = addDataListener(obj,funHandle)
            l = addlistener(obj,'dataUpdated',funHandle);
        end
    end
    
end

