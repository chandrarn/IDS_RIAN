function [tout] = imp_time_shift(tin, shot, array, ...
    probe, dir)

% this function shifts the time base of the probe signals due to time
% base differences between digitizers

shift612 = 5e-6;
shift2412 = 2.5e-6;
tout = tin;
% this will be overwritten if necessary below

if shot < 117860
    tout = tin;
    
elseif shot >= 117860 && shot <= 118389
    if strcmp(array, 'M') == 1
        if strcmp(probe, '06') == 1
            if strcmp(dir, 'R') == 1
                tout = tin - shift2412;
            end
        end
    elseif strcmp(array, 'B') == 1
        if strcmp(probe, '02') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '03') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '04') == 1
            tout = tin - shift2412;
        elseif strcmp(probe, '05') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '06') == 1
            tout = tin - shift612;
        end
    end

%     shifted probes:            
%     2412:
%     bbot probe 4 all dir
%     bmid probe 6 rad
%     612:
%     bbot probes 2,3,5,6 all dir
    
elseif shot > 118389 && shot <= 121973
    if strcmp(array, 'M') == 1
        if strcmp(probe, '08') == 1
            if strcmp(dir, 'R') == 1
                tout = tin - shift2412;
            end
        elseif strcmp(probe, '10') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '12') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '14') == 1
            tout = tin - shift612;
        elseif strcmp(probe, '17') == 1
            tout = tin - shift612;
        end
    elseif strcmp(array, 'B') == 1
        if strcmp(probe, '08') == 1
            tout = tin - shift2412;
        end
    end

%     shifted probes:            
%     2412:
%     bbot probe 8 all dir
%     bmid probe 8 rad
%     612:
%     bmid probes 10,12,14,17 all dir
    
elseif shot > 121973 && shot <= 127542
    if strcmp(array, 'M') == 1
        if strcmp(probe, '06') == 1
            if strcmp(dir, 'P') == 1
                tout = tin - shift2412;
            elseif strcmp(dir, 'T') == 1
                tout = tin - shift2412;
            end
        elseif strcmp(probe, '15') == 1
            if strcmp(dir, 'P') == 1
                tout = tin - shift2412;
            elseif strcmp(dir, 'T') == 1
                tout = tin - shift2412;
            end
        end
        if strcmp(dir, 'R') == 1
            tout = tin - shift612;
        end
    end

%     shifted probes:            
%     2412:
%     bmid probes 6,15 pol and tor
%     612:
%     bmid all rad
   

% changed the digitization rate before these shots so the time shift is
% smaller (only the 612's are off by as much as a usec)

elseif shot > 127542
    
    tout = tin;

end








