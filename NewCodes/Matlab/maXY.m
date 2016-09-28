function [ x0,y0] = maXY(Z)
% finds reasonably accurate maxima for gauss fitting ( y is dim 1)
x0 = size(Z,2)/2;
y0 = size(Z,1)/2;
% try
        if isempty(find(isnan(Z))) % skip if its a nan, return middle
            [~,y1] = max(max(Z,[],2));
            %y1 = I-1 + yBound(1); 
            [~,x1] = max(max(Z,[],1));
            %x1 = I-1 + xBound(1);
            
            if y1 ~= 1 && y1 ~= size(Z,1)
                mag = sum(Z(y1-1:y1+1,x1));
                y0 =  ( (y1-1:y1+1) * Z(y1-1:y1+1,x1) )/mag;
            end
            
            if x1 ~= 1 && x1 ~= size(Z,2)
                mag = sum(Z(y1,x1-1:x1+1));
                x0 =  ( (x1-1:x1+1) * Z(y1,x1-1:x1+1)' )/mag;
            end
%             assignin('base','Z','Z');
        else
            x0 = size(Z,2)/2;
            y0 = size(Z,1)/2;
            display('WASNAN');
        end
% catch
%     assignin('base','Z','Z');
%     display('BROKEN');
% end
%             
end