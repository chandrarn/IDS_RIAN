function [data,time] = loadCalibMovies(shot3,shot4,tOffset,trimData)

addpath('T:\PhantomMovies');
addpath('C:\Users\Rian\Dropbox\IDS_IDK');

try
     data = importdata(['Shot ' int2str(shot3) '.mat']);
     time3 = data.TimeVector;
     data = data.CineArray;
     % This is necessary for shitty, shitty back compatibility
     %data = shiftdim(data,2); % unnecessary for 161017xxx
     if ~isempty(tOffset) % apply temporal offset movie one
        data = data(tOffset(1):end,:,:); 
        time3 = time3(tOffset(1):end);
     end
     if isempty(trimData)
         trimData1 = 1:size(data,1);
     else
         trimData1=trimData;
     end
%      Unfortunately, the following all must be within the TRY because the
%      oldformat is moronic, and put time as the first column.
% The sizes here sometimes need to be manually changed, IE: avi converter
% may do the indexing differently. note 'shiftdim'
    
     [ n_time,n_pix, n_chan,] = size(data);
     data3 = zeros(length(trimData1), n_pix, n_chan);
     for i = 1:length(trimData1)
         data3(i, 1:n_pix, 1:n_chan) = squeeze(data(trimData1(i),:,:));
     end
%      data3 = shiftdim(data,3);
     data3 = data3(:, end:-1:1, :);%(trimData,:,:);
     time3=time3(trimData1);
catch
    
    data3 = importdata(['shot' int2str(shot3) '.mat']); % [counts] (time x wavelength space x channel space)
    time3 = importdata(['t' int2str(shot3) '.mat']); % [ms]
    
    if isempty(trimData)
        trimData = 1:size(data3, 1);
    end
    data3 = data3(trimData, end:-1:1, end:-1:1);
    time3 = time3(trimData);
    clear h1;
end
[n_time, n_pix, n_chan] = size(data3);
assignin('base','data3',data3);
% THE EXPLICIT TIME TRIMMING IS BECAUSE YOU DONT WANT TO FIT TO THE SECOND
% LINE

if shot4 ~= 0 % using both fibers
     try
         data = importdata(['Shot ' int2str(shot4) '.mat']);
         time4 = data.TimeVector;
         data = data.CineArray;
         
         if ~isempty(tOffset) % apply temporal offset movie one
            data = data(tOffset(2):end,:,:); 
            time4 = time4(tOffset(1):end);
         end
         %data = shiftdim(data,2); % unnecessary for 161017xxx
         if isempty(trimData)
             trimData2 = 1:size(data,1);
         else
             trimData2=trimData;
         end
%          Unfortunately, the following all must be within the TRY because the
%          oldformat is moronic, and put time as the first column.
         [n_time,n_pix, n_chan] = size(data);
         data4 = zeros(length(trimData2), n_pix, n_chan);
         for i = 1:length(trimData2)
             data4(i, 1:n_pix, 1:n_chan) = squeeze(data( trimData2(i),:,:));
         end
         time4=time4(trimData2);
%          data4 = shiftdim(data,3);
         data4 = data4(:, end:-1:1, :);%(trimData,:,:);
    catch

        data4 = importdata(['shot' int2str(shot4) '.mat']); % [counts] (time x wavelength space x channel space)
        time4 = importdata(['t' int2str(shot4) '.mat']); % [ms]

        if isempty(trimData)
            trimData = 1:size(data4, 1);
        end
        data4 = data4(trimData, end:-1:1, end:-1:1);
        time4 = time4(trimData);
        clear h1;
    end
    % Add movies together

    if length(time3) > length(time4)
                size(time3)
%         size(time4)
%         size(data3)
%         size(data4)
        data = data4 + data3(1:length(time4), :, :);
        time = time4;
        clear data3 data4 time3 time4;
    else
%         size(time3)
%         size(time4)
%         size(data3)
%         size(data4)
        data = data3 + data4(1:length(time3), :, :);
        time = time3;
        clear data3 data4 time3 time4;
    end
else % only using one fiber
    data = data3;
    time = time3(trimData);
    clear data3 time3;
end

%cd('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');

data = cast(data, 'double');
assignin('base','data',data);

end