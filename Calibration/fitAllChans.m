function peaks = fitAllChans(data, chanNums, firstCenter, lastCenter, brightWing,xWing, approxFit, force, remove)

plotFit = 1;

addpath('T:\RChandra\A-A-Ron Code\General Matlab\extrema');
addpath('/media/alfventemp/IDS/General Matlab/extrema');

% bin data into line

[n_pix, n_chan] = size(data);

lamBound = round(firstCenter(2)) - brightWing : round(lastCenter(2)) + brightWing;
chanBound = round(firstCenter(1)) - xWing : round(lastCenter(1)) + xWing;

% Catches in case brightWing extends beyond CCD

if chanBound(1) < 1
    chanBound = chanBound(find(chanBound == 1)) : chanBound(end);
end
if chanBound(end) > n_chan
    chanBound = chanBound(1):chanBound(find(chanBound == n_chan));
end

% Sum data in wavelength direction

data = sum(data(lamBound, chanBound), 1);

peaks(:, 1) = chanNums; % put channel numbers into peaks array

% Use simple interpolation to guess at peak locations in wavelength space

allChans = linspace(firstCenter(2), lastCenter(2),length(chanNums));%chanNums(end) - chanNums(1) + 1);
%size(allChans(chanNums))

%peaks(:, 3) = allChans(chanNums); % remove dead channels
peaks(:, 3) = allChans;

% Find local maxima using 'extrema' function
[xmax, imax, xmin, imin] = extrema(data); % xmax is local maxima, imax is indices
if ~isempty(remove)
    remove = remove - chanBound(1) + 1;
    chanBound
    for n = 1:length(remove)
        imax
        ind = find(remove(n) == imax);
        ind
        imax = imax([1:ind-1, ind+1:end]);
        find(imax ~= ind)
        imax
    end
end
if ~isempty(force)
    for n = 1:length(force) % index of 'data' does not start at 1. Make 'force' correspond to index of data
        force(n) = find(force(n) == chanBound);
    end
    imax = [imax, force]; % add on forced index of maxima
end
% size(imax)
% %TEST: remove all max values over number of chans
% if max(size(imax)) >max(chanNums);
%     imax=imax(1:max(chanNums));
% end
assignin('base','imax',imax);

centers = sort(imax(1:length(chanNums))); % sort indices in ascending order
amps = data(centers); % arrange maxima values to match indices
centers = chanBound(1) + centers - 1;
if centers(1) < firstCenter(1) - 1; % catch if it finds a max at the very beginning
    centers = centers(2:end);
    amps = amps(2:end);
end
if centers(end) > lastCenter(1) + 1; % catch if it finds a max at the very beginning
    centers = centers(1:end-1);
    amps = amps(1:end-1);
end
disp(['Number of extrema found: ' num2str(size(centers, 2))]);
disp(['Requested: ' num2str(length(chanNums))]);

assignin('base','centers',centers)
assignin('base','amps',amps)

if length(centers)>length(chanNums)
    centers=centers(1:length(chanNums))';
    amps=amps(1:length(chanNums));
end

% Plot summed Data

if plotFit
    S = get(0,'ScreenSize');
    fntsz = 20;

    h1 = figure('Visible','on','Name','Approximate Fit','Position',...
        [S(3)/12, S(4)/6, 5*S(3)/6 2*S(4)/3], 'Color', [1 1 1]);
    h2 = axes('Parent', h1, 'Position', [.1 .2 .8 .6], 'FontSize', fntsz);
    h3 = plot(chanBound, data, '-b');
    hold on;
    h4 = plot(centers, amps, '+r');
    for n = 1:length(centers)
        text(centers(n), 1.15*amps(n), int2str(chanNums(n)));
    end
    grid on;
    set(gca, 'XLim', [chanBound(1) chanBound(end)]);
    xlabel('Pixel Number (Real Space)');
    ylabel('Pixel Number (Wavelength Space)');
    
    % SAVE
    if approxFit.save
        fig_save = getframe(h1);
        [Xfig, mapfig] = frame2im(fig_save);
        imwrite(Xfig, [approxFit.file '.png']);
    end
end
size(peaks(:,2))
size(centers)
peaks(:, 2) = centers; % store approximate centers in peaks array

end