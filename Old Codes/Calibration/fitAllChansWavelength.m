function peaks = fitAllChansWavelength(data)

%plotFit = 1;

addpath('T:\RChandra\A-A-Ron Code\General Matlab\extrema');

% bin data into line
%% I DONT NECCESSARILY NEED TO DO THIS ( channel bounds )
%{
% [n_pix, n_chan] = size(data);
% 
% lamBound = round(firstCenter(2)) - brightWing : round(lastCenter(2)) + brightWing;
% chanBound = round(firstCenter(1)) - brightWing : round(lastCenter(1)) + brightWing;
%}
%% OR THIS ( catch chan bound exceed screem dim )
%{
% % 
% % % Catches in case brightWing extends beyond CCD
% % 
% % if chanBound(1) < 1
% %     chanBound = chanBound(find(chanBound == 1)) : chanBound(end);
% % end
% % if chanBound(end) > n_chan
% %     chanBound = chanBound(1):chanBound(find(chanBound == n_chan));
% % end
%}
%%
% Sum data in wavelength direction

%data = sum(data(lamBound, chanBound), 1);
data = sum(data, 2); % sum in the real space direction ( rows )

%peaks(:, 1) = chanNums; % put channel numbers into peaks array MAY be
%unneccessary, not sure. Depened is we want to guess where the channels are
%or not.
%% GUESS PEAK LOCATIONS, MAY NOT NEED
%{
% Use simple interpolation to guess at peak locations in wavelength space
% allChans = linspace(firstCenter(2), lastCenter(2), chanNums(end) - chanNums(1) + 1);
% peaks(:, 3) = allChans(chanNums); % remove dead channels
%}
% Find local maxima using 'extrema' function

[xmax, imax, xmin, imin] = extrema(data); % xmax is local maxima, imax is indices

peaks(:,1)=[1;2]; %first and second peak
peaks(:,2)=sort(imax(1:2)'); % wavelength coordinate of first/second peak,
%this is sorted so that we know which is the top and which is the bottom.

%% Dont need the force add, or catch for peak at begnining or end
%{
% if ~isempty(force)
%     for n = 1:length(force) % index of 'data' does not start at 1. Make 'force' correspond to index of data
%         force(n) = find(force(n) == chanBound);
%     end
%     imax = [imax, force]; % add on forced index of maxima
% end
% centers = sort(imax); % sort indices in ascending order
% amps = data(centers); % arrange maxima values to match indices
% centers = chanBound(1) + centers - 1;
% if centers(1) < firstCenter(1) - 1; % catch if it finds a max at the very beginning
%     centers = centers(2:end);
%     amps = amps(2:end);
% end
% if centers(end) > lastCenter(1) + 1; % catch if it finds a max at the very beginning
%     centers = centers(1:end-1);
%     amps = amps(1:end-1);
% end
%}

disp(['Number of Wavelength extrema found: ' num2str(size(peaks(:,2)))]);
%disp(['Requested: ' num2str(length(chanNums))]);

%% DEFINITELY DONT NEED 
%{
% Plot summed Data
% 
% if plotFit
%     S = get(0,'ScreenSize');
%     fntsz = 20;
% 
%     h1 = figure('Visible','on','Name','Approximate Fit','Position',...
%         [S(3)/12, S(4)/6, 5*S(3)/6 2*S(4)/3], 'Color', [1 1 1]);
%     h2 = axes('Parent', h1, 'Position', [.1 .2 .8 .6], 'FontSize', fntsz);
%     h3 = plot(chanBound, data, '-b');
%     hold on;
%     h4 = plot(centers, amps, '+r');
%     for n = 1:length(centers)
%         text(centers(n), 1.15*amps(n), int2str(chanNums(n)));
%     end
%     grid on;
%     set(gca, 'XLim', [chanBound(1) chanBound(end)]);
%     xlabel('Pixel Number (Real Space)');
%     ylabel('Pixel Number (Wavelength Space)');
%     
%     % SAVE
%     if approxFit.save
%         fig_save = getframe(h1);
%         [Xfig, mapfig] = frame2im(fig_save);
%         imwrite(Xfig, [approxFit.file '.png']);
%     end
% end
% peaks(:, 2) = centers; % store approximate centers in peaks array
% 
% end

%}