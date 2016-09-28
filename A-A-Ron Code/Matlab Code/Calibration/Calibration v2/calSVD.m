function[data, X, Y, n_pix, n_chan] = calSVD(shot3, shot4)



 try
     data=importdata(['Shot ' int2str(shot3) '.mat']);
     data = data.CineArray;
     clear CineArray TimeVector;
     display('Using MATLAB Cine Converter');
     %oldformat is moronic, and put time as the first column.
     [ n_pix, n_chan, n_time] = size(data);
     BDdat = NaN*zeros(n_chan * n_pix, n_time);
     for n = 1:n_time
         BDdat(:, n) = reshape(squeeze(data(:, :,n)), n_chan * n_pix, 1);
     end
 catch
    display('Using Python Cine Converter');
    %CANT LOAD MY VERSION, THE DATA IS CRAP
%     try
%         cd('T:\PhantomMovies');
%         data = importdata(['shot' int2str(shot3) '.mat']); % [counts] (time x wavelength space x channel space)
%     catch
%         cd('G:\Cines'); %If it cant load from the temp drive,
%            %copy the file over to the local drive, and read it from there.
%         data = importdata(['shot' int2str(shot3) '.mat']);
%     end
    %cd('G:\Cines'); %If it cant load from the temp drive,
    %addpath('E:\Cines\SillyFolder');
    addpath('T:\PhantomMovies');
    %copy the file over to the local drive, and read it from there.     data = importdata(['shot' int2str(shot3) '.mat']);
    % Unflip images to correct for flip during Python conversion
    data = importdata(['shot' int2str(shot3) '.mat']);
    data = data(:, end:-1:1, end:-1:1);
    [n_time, n_pix, n_chan] = size(data);
    
    BDdat = NaN*zeros(n_chan * n_pix, n_time);
    for n = 1:n_time
        BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
    end
end

% Arrange in vertical columns

%  DOESNT WORK WITH NEW CINE ARRAY FORMAT, PUT IN CATH:OLD FOMAT
% BDdat = NaN*zeros(n_chan * n_pix, n_time);
% 
% for n = 1:n_time
%     BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
% end
% clear data;

BDdat = BDdat - min(BDdat(:)); % subtract minimum value in data

% Perform SVD -----------------------

[U, S, V] = svd(BDdat, 'econ');
clear S V;

% Pull out first topo

topo1 = reshape(U(:, 1), n_pix, n_chan);
clear U;

% Make sure data is positive

if min(topo1(:)) < 0
    topo1 = -topo1;
end

[X, Y] = meshgrid(1:n_chan, 1:n_pix);

if shot4 ~= 0 % load the second shot, if it exists
    
    
%     try
%         clear data
%         cd('T:\PhantomMovies');
%         data=importdata(['Shot ' int2str(shot4) '.mat']);
%         data = data.CineArray;
%         clear CineArray TimeVector;
%       
%         %Also needs to go here because someone at NSTX is bad at life.
%         [n_pix, n_chan,n_time] = size(data);
%         BDdat = NaN*zeros(n_chan * n_pix, n_time);
% 
%         for n = 1:n_time
%             BDdat(:, n) = reshape(squeeze(data( :, :,n)), n_chan * n_pix, 1);
%         end
%     catch
%CANT USE NEW FORMAT DATA, IT DOESNT WORK. 
% move the shot files away from temp drive, matlab doesnt like it.
%          try
%              cd('T:\PhantomMovies');
%             data = importdata(['shot' int2str(shot3) '.mat']); % [counts] (time x wavelength space x channel space)
%         catch
            %%cd('G:\Cines');
            data = importdata(['shot' int2str(shot4) '.mat']);
%         end

        % Unflip images to correct for flip during Python conversion
        data = data(:, end:-1:1, end:-1:1);
        [n_time, n_pix, n_chan] = size(data);
        BDdat = NaN*zeros(n_chan * n_pix, n_time);

        for n = 1:n_time
            BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
        end
%     end
    %data = importdata(['shot' int2str(shot4) '.mat']); % [counts] (time x wavelength space x channel space)
    %addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Calibration\Calibration v2');
    addpath('/media/alfventemp/RChandra/A-A-Ron Code/Matlab Code/Calibration/Calibration v2');
    addpath('T:\IDS\Calibration');
    % Unflip images to correct for flip during Python conversion

    %data = data(:, end:-1:1, end:-1:1);

    % Arrange in vertical columns
% 
%     [n_time, n_pix, n_chan] = size(data);
%     BDdat = NaN*zeros(n_chan * n_pix, n_time);
% 
%     for n = 1:n_time
%         BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
%     end
    clear data;

    BDdat = BDdat - min(BDdat(:)); % subtract minimum value in data

    % Perform SVD -----------------------

    [U, S, V] = svd(BDdat, 'econ');
    clear S V;

    % Pull out first topo

    topo2 = reshape(U(:, 1), n_pix, n_chan);
    clear U;
    
    % Make sure data is positive

    if min(topo2(:)) < 0
        topo2 = -topo2;
    end

    % Add topos from both shots together
    size(topo1)
    size(topo2)
    data = topo1 + topo2;
    clear topo1 topo2;
else
    addpath('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');
    data = topo1;
    clear topo1;
end






