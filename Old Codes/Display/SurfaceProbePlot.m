% Display surface probe data
% One page per probe number, 1-10, for large and small cones
% display torroidal and poloidal in columns, of decreasing frequency.

%% NOTE: THERE IS AN ISSUE WITH XLABEL ON THE LAST PLOT: NO IDEA WHY

cd('T:\IDS\Display\Filters'); % location of surface probe finder, a la brian victor

clear all; close all;

[bsp, phi, theta, node_save] = sp_b_in_updated(-.0005:.000002:.004); % get the data, over some time base

names= fieldnames(bsp); % get the names of the probes

isSmall = strfind(names,'S'); % finds the ones that have S. L goes to []
isSmall = ~cellfun('isempty',isSmall); % This is magic. Boolean array of Small cone probes
Small = names(isSmall);
Large = names(~isSmall);
% Temp = 0; % Convert to string array
% for i = 1:length(Small)
%     Temp(i) = Small{i};
% end
% Small = Temp; Temp = 0;
% for i = 1:length(Large)
%     Temp(i) = Large{i};
% end
% Large = Temp;

%Initialize the probe struct
SortedNames = cell(length(names),1);
SortedNames(1:length(Small)) = Small;
SortedNames(length(Small)+1:end) = Large;

Probes.Small = [0];
Probes.Large = [0];
ItrSmall = 1;% initialize iterators
ItrLarge = 1;
itr = 1;
deadCounter = 1;

%figure stuff
S=get(0,'ScreenSize');
analysisHeight = S(4)-110; % all the sky
figureWidth = (S(3)-12)/2; % half the sky
r = figure('Visible','on','Name','Surface Probes','Position',...
    [5 35 figureWidth analysisHeight],'Color',[1 1 1],'Units','Normalized');
for i = 1:length(SortedNames) % figure out where to put all the probes 
    
    h(itr) = subplot(10,2,itr);
    
    %set(r,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
    a(i,:) = get(h(itr),'Position');
    if itr ~= 1 % because subplot is a moron.
        set(h(itr),'Position',[ a(i,1:3),0.0603703703703704]);
    end
    eval(['hPlot(i)=plot(h(itr),[-.5:.002:4],bsp.' SortedNames{i} '.*100);']);
    
    %store nans
    if eval(['any(isnan(bsp.' SortedNames{i} '))']) % if dead probe
        deadProbes(deadCounter,1:2) = a(i,1:2);
        deadCounter=deadCounter+1;
    end
    grid on;
    %legend(char({[name(5:7) ', '] , [name(9:end) char(176)]}));
%     temp(i,:) = [a(1)+.275,a(2)+.005,.05,.05];
%      hold on;
%     an(i)=annotation('textbox',[a(1)+.275,a(2)+.005,.05,.05],'String', ...
%         {[name(5:7) ', '] , [name(9:end) char(176)]},'BackgroundColor',[1 1 1 ]);
     %inputdlg('TEST');
    hold on;
    set(h(itr),'XLim',[-.5,3.5]);
    
    %Special plots:
    if itr == 19 | i == length(names)-1 
        xlabel('Time [ms]');
        ylabel('B [mT]');
        
    elseif itr ==20 | i == length(names)
        % labels
        xlabel('Time [ms]');
        ylabel('B [mT]');
        
        % stupid friggen annotations loop
        for j = i-itr+1:i
            name = SortedNames{j};% CURLY BRACKETS IMPORTANT
            annotation('textbox',[a(j,1)+.275,a(j,2)+.005,.05,.05],'String', ...
        {name(5:7) , [name(9:end) char(176)]},'BackgroundColor',[1 1 1 ]);
        end
        
        % titles
        annotation('textbox',[.25,.93,.1,.025],'String', ...
            {'   Poloidal'},'BackgroundColor',[1 1 1 ]);
        annotation('textbox',[.7,.93,.1,.025],'String', ...
            {'   Toroidal'},'BackgroundColor',[1 1 1 ]);
        %inputdlg('TEST');
        
        %label dead probes
        if deadCounter ~= 1 % if deadprobes exist
        for k = 1:size(deadProbes,1)
            annotation('textbox',[deadProbes(k,1)+.1-.025,deadProbes(k,2)+.02,.15,.025],...
                'String',{'   DEAD PROBE'},'BackgroundColor',[1 1 1 ]);
        end
        end
        clear deadProbes;
        deadCounter = 1;
        
        itr=0;
        saveas(r,['T:\IDS\Analysis Repository\SurfacePlot129817' num2str(i)],'eps'); %HAS TROUBLE WITH FORMATTING
%         fig_save = getframe(r);
%         [Xfig, mapfig] = frame2im(fig_save);
%         imwrite(Xfig, ['T:\IDS\Analysis Repository\SurfacePlot' num2str(i) '.bmp']);
        pause(1);
        %close all;
        r = figure('Visible','on','Name','Surface Probes','Position',...
    [5 35 figureWidth analysisHeight],'Color',[1 1 1]);
    else
        %set(gca,'XTick',[]);
        set(gca,'XTickLabel',[]);
    end
    itr=itr+1;
end
