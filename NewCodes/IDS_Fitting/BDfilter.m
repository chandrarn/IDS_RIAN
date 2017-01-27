function [data, time] = BDfilter(shot, nBDmodes, timeBound, param, s)
% Aaron Hossack
%
% The '3' version does a better job- now, tight bounds are set around the
% toroidal and poloidal data, they are individually BD'd, and then placed
% onto an array of zeros.  Thus, the movie is the same size as the original
% so the channel numbers match.

plotFits = 0;
BDwing = param.yWing ;%+ 5  % number of points to either side of center in wavelength space to include in BD domain

modes = 1:nBDmodes; % array of mode numbers to recombine and add together

first = param.peaks(1, 2);
last = param.peaks(find(param.peaks(:, 1) <= 36, 1, 'last'), 2);
chan_bound_t = (floor(first)):(ceil(last)+10);

first = param.peaks(find(param.peaks(:, 1) >= 37, 1), 2);
last = param.peaks(end, 2);
chan_bound_p = (floor(first)):(ceil(last)+0);%+10

%remove -3 if you see it 23/2/16
%lam_bound = round(param.Center(1, 1)) - BDwing : round(param.Center(1, end)) + BDwing +18
lam_bound = round(param.Center(1, 2))-BDwing -28  : round(param.Center(1, 2))+BDwing - 28 ;
% Since 'Center' is ordered longes to shortest the indexing is a little
% weird.

%%
% Load Data Array

addpath('T:\PhantomMovies');
%Try to import the files from the new format, if possible
try
    dummy = importdata(['Shot ' int2str(shot) '.mat']);
    data = dummy.CineArray;
    time = dummy.TimeVector;
    data = data(end:-1:1, :,:);%precut time
    newCineType = 1;
    clear dummy;
catch

    if ~isempty(s.sim)
        shot = str2num([int2str(s.sim) int2str(shot)]); % FUCK YEAH
    end

        data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
        time = importdata(['t' int2str(shot) '.mat']); % [ms]
        data = data(:, end:-1:1, end:-1:1); % removed x and y flip
        newCineType=0;
end
assignin('base','time',time)
%cd('T:\RChandra\A-A-Ron Code\Matlab Code\Core Fitting Codes\BDfilter');
[n_time_o, n_pix_o, n_chan_o] = size(data);


% Account for the possibility of no wavelength bounds
if isempty(lam_bound)
    lam_bound = 1:n_pix_o;
end

% Account for empty channel bound
if isempty(chan_bound_t)
    chan_bound_t = 1:n_chan_o;
end
if isempty(chan_bound_p)
    chan_bound_p = 1:n_chan_o;
end

% BD specific time bound
if isempty(timeBound)
   timeBound = 1:length(time)-1;
end
size(time)
size(timeBound)
time = time(timeBound);



for m = 1:2
    if m == 1 % toroidal fiber
        chan_bound = chan_bound_t;
    else % poloidal fiber
        chan_bound = chan_bound_p;
    end
    
   
%     % BD specific time bound
%     if isempty(timeBound)
%         timeBound = 1:length(time);
%     end
%     time = time(timeBound);

    %isolate time in new/old filetypes
    if(~newCineType)
        
        % Record size for recombining
        [n_time_o, n_pix_o, n_chan_o] = size(data)
        timeBound
        lam_bound
        chan_bound
        data2 = data(timeBound, lam_bound, chan_bound);
        
        % Arrange in vertical columns
        [n_time, n_pix, n_chan] = size(data2);
        BDdat = NaN*zeros(n_chan * n_pix, n_time);
        for n = 1:n_time
            BDdat(:, n) = reshape(squeeze(data2(n, :, :)), n_chan * n_pix, 1);
        end
    else
        
        [ n_pix_o, n_chan_o, n_time_o] = size(data) % DONT MODEFY DATA IN THE LOOP:
        % BECAUSE IT GOES TWICE, IT WILL TRY TO TRIM IT AGAIN, WHICH IS BAD
        data2=data(lam_bound,chan_bound,timeBound);
        % Arrange in vertical columns
        [n_pix, n_chan, n_time] = size(data2);
        BDdat = NaN * zeros(n_chan * n_pix, n_time);
        for n = 1:n_time
            BDdat(:, n) = reshape(squeeze(data2(:, :, n)), n_chan * n_pix, 1);
        end
    end

    
    
    % Subtract off minimum value to reduce noise floor
    minVal = min(BDdat(:));
    
    if plotFits
        disp(['Min = ' num2str(minVal)]);
        disp(['Median = ' num2str(median(BDdat(:)))]);
        disp(['Mean = ' num2str(mean(BDdat(:)))]);
        disp(['Std = ' num2str(std(BDdat(:)))]);
    end
    
    BDdat = BDdat - minVal; % subtract minimum value in data
    
    %% Perform SVD -----------------------

    [U, S, V] = svd(BDdat, 'econ');
    Ak = diag(S);
    
    % Rearrange data as it was before
    %% NOTE: LEAVING TOPOS IN (T,X,Y) FORMAT FOR CONSISTANCY.I STILL THINK THAT THIS
    % FORMAT IS ASCININE, BUT WHATEVER
    topos = zeros(n_time, n_pix, n_chan);
    
    for n = 1:n_time % actually a loop over topo modes
        topos(n, :, :) = reshape(U(:, n), 1, n_pix, n_chan);
    end
    
    %% Plotting
    if plotFits
        [X, Y] = meshgrid(chan_bound, lam_bound);
        
        S1 = get(0,'ScreenSize');
        lnwdth = 2;
        fntsz = 20;
        
        h1 = figure('Visible','on','Name','Topo 1','Position',...
            [S1(3)/4, S1(4)/6, S1(3)/2 2*S1(4)/3], 'Color', [1 1 1]);
        h2 = axes('Parent', h1, 'Position', [.1 .1 .8 .8], 'FontSize', fntsz);
        h3 = surf(X, Y, Ak(1)*squeeze(topos(1, :, :)));
        hold on;
        shading interp;
        colormap jet;
        colorbar;
        grid on;
        view([0 90]);
        title(['First Topo, shot ' num2str(shot)]);
        xlabel('Pixel Number (real space)');
        ylabel('Pixel Number (wavelength space)');
        set(h3, 'LineWidth', lnwdth);
        set(gca, 'FontSize', fntsz, 'XLim', [1 n_chan], 'YLim', [1 n_pix]);
    end
    
    %% Save Data for Recombining
    if m == 1
        Ak_t = Ak(modes); % save weights
        figure; plot(Ak_t,'-*','linewidth',2);title('Toroidal BD modes')
        V_t = V(:, modes); % save chronos
        topos_t = topos(modes, :, :); % save topos
        assignin('base','Ak_t',Ak);
    else
        Ak_p = Ak(modes); % save weights
        V_p = V(:, modes); % save chronos
        topos_p = topos(modes, :, :); % save topos
        assignin('base','Ak_p',Ak);
    end
    clear Ak V topos
    
end
%% Recombine Selected Modes into False Shot ASSUME RECOMBINE IS 1


assignin('base','data2',data2)
% NOTE: CHANGED N_TIME_O TO N_TIME: THIS SHOULD ONLY MATTER WHEN TIMEBOUND IS ACTIVE
data = zeros(n_time, n_pix_o, n_chan_o); % THIS SHOULD BE COMPATABLE WITH NEW CINE ARQUETECTURE
data2 = data; % initialize both to zeros

for n = 1:length(modes)
    for m = 1:size(V_t, 1) % loop over time
        data2(m, lam_bound, chan_bound_t) = Ak_t(n) * V_t(m, n) * topos_t(n, :, :);
        data2(m, lam_bound, chan_bound_p) = Ak_p(n) * V_p(m, n) * topos_p(n, :, :);
    end
    data = data + data2;
end

%toc
%clear data2 Ak_t V_t topos_t Ak_p V_tp topos_p;
%assignin('base','Ak_t',Ak_t);assignin('base','Ak_p',Ak_p);
figure; 
surf(squeeze(sum(data, 1)./size(data,1))); 
shading interp;
view([ 0 90]);



