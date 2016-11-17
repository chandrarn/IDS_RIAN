% Analysis of Modes and phases from coherant fiber movie

shot = 160412022; % Deuterium
load(['T:\PhantomMovies\Shot ' num2str(shot) '.mat']);

% movie:
v=VideoWriter('Coherant_Plasma_Movie2.mp4','MPEG-4');
%v.LosslessCompression=true;
v.Quality=100;
v.FrameRate=2;
open(v);
h=figure;
for i = 500:600
    i
    surf(squeeze(CineArray(:,:,i))); shading interp; view([ 0 90]);
    title(['Coherant Movie Time Point: ' num2str(TimeVector(i).*1e-3) '[ms]']);
    colorbar; set(gca,'clim',[-50,3000]);
    pause(.75);
    writeVideo(v,getframe(h));
end
close(v);


redColorMap = [linspace(1, 0, 124), linspace(0,1, 132)];
colorMap = [redColorMap; zeros(1, 256); zeros(1, 256)]';

% A injector
for i = 1:120
for j = 1:113
signal=squeeze(CineArray(i,0+j,500:600));
time=TimeVector(500:600).*1e-6;
guess=[mean(CineArray(i,0+j,500:600)),std(CineArray(i,0+j,500:600)),0,14500];
[param(i,j,1:4)]=sine_fit(time',signal',[nan,nan,nan,14500],guess,0);
end
i
end
figure; surf(mod((squeeze(param(:,:,3)-pi).*(param(:,:,2)<0)) + (squeeze(param(:,:,3)).*(param(:,:,2)>0)),2*pi));
shading interp; view([ 0 90]);colormap(colorMap);colorbar;
title('Phase of A Injector Mouth');

% Modes
modes=10;
[n_pix, n_chan, n_time] = size(CineArray(:,1:112,:));
BDdat = NaN * zeros(n_chan * n_pix, n_time);
for n = 1:n_time
    BDdat(:, n) = reshape(squeeze(CineArray(:, 1:112, n)), n_chan * n_pix, 1);
end
 minVal = min(BDdat(:));
    
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
data2 = zeros(n_time, n_pix, n_chan);
for n = 1:length(modes)
    for m = 1:size(V, 1) % loop over time
        data2(m, :, :) = Ak(n) * V(m, n) * topos(n, :, :);
    end
end
% Plot first mode
figure;
surf(squeeze(Ak(1) * V(i, 1) * topos(1, :, :))); shading interp; view([ 0 90]);
colormap(parula); title('First BD Mode, A Inj Mouth')




% Full volume
for i = 1:120
for j = 1:113
signal=squeeze(CineArray(i,111+j,500:600));
time=TimeVector(500:600).*1e-6;
guess=[mean(CineArray(i,111+j,500:600)),std(CineArray(i,111+j,500:600)),0,14500];
[param2(i,j,1:4)]=sine_fit(time',signal',[nan,nan,nan,14500],guess,0);
end
i
end
% Plot Phase
figure; surf(mod((squeeze(param2(:,:,3)-pi).*(param2(:,:,2)<0)) + (squeeze(param2(:,:,3)).*(param2(:,:,2)>0)),2*pi));
shading interp; view([ 0 90]);colormap(colorMap)
colorbar;
title('Phase of Confinement Volume');
% Modes
modes=10;
[n_pix, n_chan, n_time] = size(CineArray(:,112:end,:));
BDdat = NaN * zeros(n_chan * n_pix, n_time);
for n = 1:n_time
    BDdat(:, n) = reshape(squeeze(CineArray(:, 112:end, n)), n_chan * n_pix, 1);
end
 minVal = min(BDdat(:));
    
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
data2 = zeros(n_time, n_pix, n_chan);
for n = 1:length(modes)
    for m = 1:size(V, 1) % loop over time
        data2(m, :, :) = Ak(n) * V(m, n) * topos(n, :, :);
    end
end
% Plot first mode
figure;
surf(squeeze(Ak(1) * V(i, 1) * topos(1, :, :))); shading interp; view([ 0 90]);
colormap(parula);title('First BD Mode, Confinemet Volume');

