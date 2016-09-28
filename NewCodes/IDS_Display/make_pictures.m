% this program reads a *.mat file and makes a movie
clear all; close all;

cd('T:\PhantomMovies');

disp_emission = 0; % plots the total signal
disp_1D = 0; % plots the 1-D spectrum (adds channels together)
disp_oneFrame = 1; % plots a frame at a particular time
disp_frameSum = 0; % display a sum of all frames in the movie
fig4LLim = 2; % frequently changed to account for sum of pixels in movie
fig4ULim = 30;
n_start = 1;
n_end = 350;

skip = 1; % number of frames to skip to keep data size reasonable (1 = none)
lamYScale = 0; % if 0, display pixel numbers on Y scale.  If set to 1, PIX_SP 
% CENTER, and WAVELENGTH must exist in the ANALYSIS tree for the shot of interest.

shot = 13052101; %input('What shot? ');
% [a, status] = mdsopen('landau.hit::hitsi', shot);

ph_cam = importdata(['shot' int2str(shot) '.mat']);
t_cam = importdata(['t' int2str(shot) '.mat']);

cd('S:\MC_IDS\Matlab Code\Display');

% [t.iyinj, dt.iyinj, sig.iyinj] = gen_data_in('\i_inj_y');
% [t.vyinj, dt.vyinj, sig.vyinj] = gen_data_in('\v_inj_y');
% [t.psiyinj, dt.psiyinj, sig.psiyinj] = gen_data_in('\psi_inj_y');

dph_cam = double(ph_cam(1:skip:end, :, :));
t_cam = t_cam(1:skip:end);

if disp_emission
    emiss = sum(sum(dph_cam, 3), 2);
    figure(1)
    plot(emiss); % plots as function of INDEX
%     plot(t_cam, emiss); % plots as function of TIME
    title('Emission over Time');
end
    

% want to subtract off differences between pixels
% this provides the option to do a zero subtraction
% it is best if this is already done with the camera,
% but sometimes it wasn't done

sub_b = 0;

if sub_b == 1;

    mean_length = 50;
    dph_zero = mean(dph_cam(2:mean_length, :, :));

    nph_cam = zeros(size(dph_cam));
    for j = 1: length(dph_cam(:,1,1))
        nph_cam(j,:,:) = dph_cam(j,:,:) - dph_zero;
    end
    % inverts pic
    inv_cam = nph_cam(:,end:-1:1,end:-1:1);

elseif sub_b ~= 1
    % this inverts the pic
    inv_cam = dph_cam(:,end:-1:1,end:-1:1);
end

clear dph_cam;

max_cam = max(max(max(inv_cam)));
nor_cam = inv_cam/max_cam;

ZI_cam = nor_cam;
clear nor_cam;

% movie start time
ts = 0;
% movie end time
te = 1.253;

cam_range = max([length(ZI_cam(1,:,1)) length(ZI_cam(1,1,:))]);

if lamYScale
    % Pull in data from tree
    mdsopen('landau.hit::hitsi', shot);
    CENTER = mdsvalue('\IDS_CENTER');
    PIX_SP = 1e9*mdsvalue('\IDS_PIX_SP');
    LAMBDA = 1e9*mdsvalue('\LINE_LAMBDA');
    mdsclose();
    
    % Select median value since tree values contain every channel
    CENTER = median(CENTER);
    PIX_SP = median(PIX_SP);
    % Select calibration Lambda
    LAMBDA = LAMBDA(1);
    PIX_SP = PIX_SP*1.5;
    ytick = [0:cam_range] .* PIX_SP + LAMBDA - CENTER * PIX_SP;
else
    ytick = [0:cam_range]; % default setting, plots by pixel
end

[aa, Is] = min((t_cam - ts).^2); % start index
[aa, Ie] = min((t_cam - te).^2); % start index

xipmin = ts - .05;
xipmax = te + .05;
yipmin = -15;
yipmax = 15;

xt1 = round(100*(xipmin + .1*(xipmax - xipmin)))/100;
xt2 = round(100*(xipmin + .5*(xipmax - xipmin)))/100;
xt3 = round(100*(xipmin + .9*(xipmax - xipmin)))/100;

scrsz = get(0, 'ScreenSize');
    
if disp_oneFrame
    h1 = figure(2);     %[Left edge, Bottom edge, Width, Height]
    set(h1, 'Position', [1 scrsz(4)/10 scrsz(3)/3 scrsz(4)/2], 'Color', [1 1 1]);
    hold('all');
    clf;

    k2 = 0;

    for k = Ie
        k2 = k2 + 1;
        axes1 = axes('Parent', h1, 'YTick', [0:cam_range], 'XTick', [0:cam_range],...
            'Position', [0.15 0.15 0.65 0.75], 'FontSize', 14,...
            'XColor', [0 0 0], 'YColor', [0 0 0], 'YTickLabel', ytick);
        box('on');
        hold('all');
        surf(squeeze(ZI_cam(k, :, :)));
        shading interp;
        caxis([0, 1]);
        xlim([0 cam_range]);
        ylim([0 cam_range]);
        colormap jet;
        colorbar('peer', axes1, [0.8068 0.5814 0.05222 0.3352],...
            'YTick', [1, .1],...
            'YTickLabel', {'Bright', 'Dark'});
        title(['Single Frame'], ...
            'FontSize', 18, 'interpreter', 'latex');
        text(.2, .05, ['time = ' num2str(t_cam(k))]);

        box('on');
        hold('all');

    end
end

%% Adding Frames
if disp_1D
    figure(3)
    plot(squeeze(sum(sum(ZI_cam, 3), 1)));
    hold on;
    title({'One Dimensional Spectrum'; '(sum over time and spatial dimensions)'});
    ylabel('counts');
    if lamYScale
        axes('XTickLabel', ytick);
    else
        xlabel('Pixel Number');
    end
end

%% End
% filename = ['S:\Fast_Camera\Movies\pos_spi_b_' num2str(shot)];
% disp(['Writing File to: ' filename])
% movie2avi(camplot, filename, 'fps', 10, ...
%     'compression', 'Cinepak', 'quality', 100); %'compression','None');

%% Sum of all frames

if disp_frameSum
    h1 = figure(4);     %[Left edge, Bottom edge, Width, Height]
    set(h1, 'Position', [1 scrsz(4)/10 scrsz(3)/3 scrsz(4)/2], 'Color', [1 1 1]);
    hold('all');
    clf;

    axes1 = axes('Parent', h1, 'YTick', [0:cam_range], 'XTick', [0:cam_range],...
        'Position', [0.15 0.15 0.65 0.75], 'FontSize', 14,...
        'XColor', [0 0 0], 'YColor', [0 0 0], 'YTickLabel', ytick);
    box('on');
    hold('all');
    surf(squeeze(sum(ZI_cam(n_start:n_end, :, :), 1)));
    shading interp;
    caxis([fig4LLim, fig4ULim]);
    xlim([0 cam_range]);
    ylim([0 cam_range]);
    colormap jet;
    colorbar('peer', axes1, [0.8068 0.5814 0.05222 0.3352],...
        'YTick', [1, .1],...
        'YTickLabel', {'Bright', 'Dark'});
    title(['Sum of Frames in Selected Time Span'], ...
        'FontSize', 18, 'interpreter', 'latex');
%     text(.2, .05, ['time = ' num2str(t_cam(k))]);
    box('on');
    hold('all');
end





