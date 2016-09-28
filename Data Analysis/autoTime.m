% auto trim time for batch correct
function [startTime, endTime] = autoTime(data,numThresh,countThresh,frameThresh)

%  loop through data, if numThresh points are greater than countThresh,
%  for frameThresh frames, then trim to that start and end point

startTime = 0;
endTime = 1;
Strig = 0; % Start at zero, if reach cutoff, start incrementing. 
Etrig = 0; % 
% if ten successive frames are above the counter, 
for n = 1:size(data,1)
    
    temp = squeeze(data(n,:,:));
    if Strig <=frameThresh % if we're looking for the start, and current frame and succeeding frame cross threshold:
        if (length(find(temp>=countThresh))>=numThresh)
            if Strig==frameThresh % if we have just reached the trigger, make this the startpoint
                startTime = n-numThresh-30;
            end
            Strig = Strig+1; % if we are crossing the threshold again, add to trigger
        else % if it falls below the threshold before it's hit ten above frames in a row
            Strig = 0; % reset
        end
    end
     % if looking for end, and curr frame and next frame fall below threshold
     if Etrig <=frameThresh && startTime > endTime
        if ~(length(find(temp>=countThresh))>=numThresh) 
            if Etrig ==frameThresh
                endTime = n -numThresh+15;
                startTime
                endTime
                output = data(startTime:endTime,:,:);
                break
            end
            Etrig = Etrig +1;
        else
            Etrig = 0;
        end
     end
end

            
    