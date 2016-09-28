% average  in space to remove Nans
% data usually more unifor in space than in time

function [data] = averageNans(data)
    for pass = 1:2 % tends to work better if we do it twice.
        if pass == 1
            tRange = 1:size(data,1);
            sRange = 2:size(data,2)-1;
        else
            tRange = size(data,1):-1:1;
            sRange = size(data,2)-1:-1:2;
        end
        for i = tRange % time
            if isnan(data(i,1))%chan 1
                data(i,1) = nanmean([data(i,2),data(i,3)]); % take the mean, exclude nans
            end
            for j = sRange % space
                if isnan(data(i,j))
                    data(i,j) = nanmean([data(i,j-1),data(i,j+1)]); % take the mean, exclude nans
                end
            end
            if isnan(data(i,end))%chan 1
                data(i,end) = nanmean([data(i,end-1),data(i,end-2)]); % take the mean, exclude nans
            end
        end
    end
    % final pass, move through time
    for i = 2:size(data,1)-1
        for j = 1: size(data,2)
            if isnan(data(i,j))
                data(i,j) = nanmean([ data(i-1,j), data(i+1,j)]);
            end
        end
    end
end
            