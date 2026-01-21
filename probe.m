classdef probe < handle
    properties(Constant)
        % ****The port for the arduino depends on the specific computer****.
         arduino = arduino('COM10','UNO'); %port for my Windows PC
        % arduino = arduino('/dev/cu.usbmodem14111','uno'); % port for my MacBook
    end 
    properties
        servobottom;
        servotop;
        position_bottom;
        position_top; 
    end
    methods
        
        % This initializes the probe object calibrated for current probe
        % characteristics. The Min and Max Pulse Duration are set to
        % provide proper scaleing for 0-180 degree rotation of servos.
        % The initializer sets the probe to its 0,5 point
        function obj = probe
            obj.servobottom = servo(obj.arduino, 'D9', 'MinPulseDuration', 0.4423e-3, 'MaxPulseDuration', 2.156e-3);
            obj.servotop = servo(obj.arduino, 'D8', 'MinPulseDuration', 0.486e-3, 'MaxPulseDuration', 2.21e-3);
            obj.position_bottom = 1.0;
            obj.position_top = 1.0;
            writePosition(obj.servobottom, obj.position_bottom);
            writePosition(obj.servotop, obj.position_top);
%             pause
%                         obj.position_bottom = 1;
%             obj.position_top = 0.5;
%             writePosition(obj.servobottom, obj.position_bottom);
%             writePosition(obj.servotop, obj.position_top);
% pause
%             for i=0:2:180
%             obj.position_bottom = i/180;
%             i
%             writePosition(obj.servobottom, obj.position_bottom);  
%             pause
%             end
%                         for i=1:2:180
%             obj.position_bottom = i/180;
%             i
%             writePosition(obj.servobottom, obj.position_bottom);  
%             pause
%             end
        end
        
        % moves object to a new location based on input. Input position
        % must be between 0 and 1, given by the matrics S1P and S2P
        % provided. The movement is done in small increments to reduce any
        % abrupt movements of the e-field probe.
        function move(obj, new_pos_bottom, new_pos_top)
            new_pos_bottom = new_pos_bottom;
            if new_pos_bottom > 1 || new_pos_bottom < 0 || new_pos_top > 1 || new_pos_top < 0
                error('Position must be between 0 and 1 for Arduino servo motors')
            else
new_pos_bottom=interp1([0 30 60 90 120 150 180]/180,[0 33 63 93 123 151 180]/180,new_pos_bottom);
new_pos_top=interp1([0 30 60 90 120 150 180]/180,[0 30 60 90 120 150 180]/180,new_pos_top);
                obj.position_bottom = new_pos_bottom;
                obj.position_top = new_pos_top;
                writePosition(obj.servobottom,obj.position_bottom);
                writePosition(obj.servotop,obj.position_top);

%     obj.position_bottom = readPosition(obj.servobottom);
%     obj.position_top = readPosition(obj.servotop);
%     
%                 writePosition(obj.servobottom,obj.position_bottom);
%                 writePosition(obj.servotop,obj.position_top);
                %                 stepb=0.005;
%                 stept=0.005;
%                 if  new_pos_bottom < obj.position_bottom
%                     stepb = stepb*-1;
%                 end
%                 if new_pos_top < obj.position_top
%                     stept = stept*-1;
%                 end
%                 for i=obj.position_bottom:stepb:new_pos_bottom
%                     writePosition(obj.servobottom,i);
%                 end
%                 for i=obj.position_top:stept:new_pos_top
%                     writePosition(obj.servotop,i);
%                 end
%                 obj.position_bottom = new_pos_bottom;
%                 obj.position_top = new_pos_top;
            end
        end 
    end
end