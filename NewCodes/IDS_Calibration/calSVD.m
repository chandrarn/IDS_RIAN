function[data, X, Y, n_pix, n_chan] = calSVD(shot1, shot2)

addpath('T:\PhantomMovies');

try
    temp=importdata(['Shot ' int2str(shot1) '.mat']); % This may fault on old matlab
    data = temp.CineArray;
    %data = data(end:-1:1, :,:); % need to flip for consistancy with Batch Correct 
    % Both could remove this, but it makes the data match visually with the
    % PCC program
    clear temp;
    [n_time, n_pix, n_chan] = size(data);
catch
    data = importdata(['shot' int2str(shot1) '.mat']); % [counts] (time x wavelength space x channel space)
    
    % Unflip images to correct for flip during Python conversion
    data = data(:, end:-1:1, end:-1:1);
end

% Arrange in vertical columns

[n_time, n_pix, n_chan] = size(data);
BDdat = NaN*zeros(n_chan * n_pix, n_time);

for n = 1:n_time
    BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
end
clear data;

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

if shot2 ~= 0 % load the second shot, if it exists

    try
        temp=importdata(['Shot ' int2str(shot2) '.mat']); % This may fault on old matlab
        data = temp.CineArray;
        %data = data(end:-1:1, :,:);
        clear temp;
    catch
        data = importdata(['shot' int2str(shot2) '.mat']); % [counts] (time x wavelength space x channel space)

        % Unflip images to correct for flip during Python conversion
        data = data(:, end:-1:1, end:-1:1);
    end

    % Arrange in vertical columns

    [n_time, n_pix, n_chan] = size(data);
    BDdat = NaN*zeros(n_chan * n_pix, n_time);

    for n = 1:n_time
        BDdat(:, n) = reshape(squeeze(data(n, :, :)), n_chan * n_pix, 1);
    end
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

    data = topo1 + topo2;
    clear topo1 topo2;
else
%     cd('S:\MC_IDS\Matlab Code\Calibration\Calibration v2');
    data = topo1;
    clear topo1;
end






