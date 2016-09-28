%function PlayMovie();

%% CREATE A MOVIE OF THE CINE AND BD FILTERED DATA, PLAY THEM

% load the array from new or old cine format
% plot frames, get frames 
% load bd filtered array
% plot frames, get frames

close all;
cd('T:\PhantomMovies');
%Try to import the files from the new format, if possible
shot=129819;
try
    dummy = importdata(['Shot ' int2str(shot) '.mat']);
    data = dummy.CineArray;
    time = dummy.TimeVector;
    data = data(end:-1:1, :,: );%flip y axis ( flip upside down )
    clear dummy;
catch
    data = importdata(['shot' int2str(shot) '.mat']); % [counts] (time x wavelength space x channel space)
    time = importdata(['t' int2str(shot) '.mat']); % [ms]
    data = data(:, end:-1:1, end:-1:1);
    data = reshape(data,size(data,2),size(data,3),[]); % WE'RE DOING IT MY WAY
end
size(data)
shot=12981910;
cd('T:\IDS\Data Repository');
eval(sprintf('load(''dat%i'');', shot));
filtered=dat.raw;
size(filtered)
%loop through the matrix, plot, get frame from plot
%BDmovie=0;
%OGmovie=0;
S=get(0,'ScreenSize');
h1 = figure('Visible','on','Name','Positive/Negative Comparison','Position',...
    [5 35 S(3)-12 S(4)-110],'Color',[1 1 1]);
for i = 1:size(data,3)
    
    surf(squeeze(dat.raw(i,:,end:-1:1))');
    shading interp;
    colormap jet;
    view([90,90]);
    set(gca, 'YLim',[0,300],'XLim',[30,45] );
    BDmovie(i)=getframe;
end

for i = 1:size(data,3)
    surf(squeeze(data(:,end:-1:1,i))');
    shading interp;
    colormap jet;
    view([90,90]);
    set(gca, 'YLim',[0,350],'XLim',[0,95] );
    OGmovie(i)=getframe;
end
movie(OGmovie);
movie(BDmovie);