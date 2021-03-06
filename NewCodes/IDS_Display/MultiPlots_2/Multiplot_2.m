%{ 
MultiPlots Improved:
Modularized the input settings, and fiber flip

The SaveDat structure is input into 
%}

clear all;
%close all;
clc;
%% Settings
% Pick set of shots and analysis settings
% n = 1;% 160728013, 160728012
% n = 7;% 129499, 129496
% n=12; %170518025;,170518028
n = 8; % 129499

[in,timebound,chan_range,xlim,supress] = in_settings(n);

% Pick plot type
plt.Type = 1; % Velocity
%plt.Type = 2; % Temperature
%plt.Type = 3; % Intensity
%plt.Type = 4; % Displacement

% Which plots to make, other settings
saving = 0; % Save the Figues
plt.FFT = 0; % Plot the FFT Power Spectrum
plt.Currents = 1; % Plot the currents at the bottom of the data
plt.SanityPhase = 0; % Plot the phase analysis stages
plt.Averages = 0; % unused
plt.compactCurrents = 1; %
plt.Tor = 0; % The line plot will show the differenve between fibers
plt.ExplainReconst=0;% will plot the data and reconstruction from line n to a seperate plot
driveDirection=1*(plt.Type==1); % upper fiber dropped by -pi, phase now references "positive toroidal flow magnitude"
flipLoImpact=[47].*((in(1).shot>229499)&&(in(1).shot<170501000)); 
% correct for lower fiber being rotated 180 about impact 30 
% the bracketed number is the pre-trimmed index of this impact. 
plt.Error = 1;
plt.includeTemp = 1; % include temperature in the velocity plot
CutPow  = .1; % FFT Power Cuttoff Percentage. 

% Spacing for line plot
velSpace = 15; % km/s
intSpace = 20; % arb.u.
tempSpace = 20; % eV
% Figure stuff
fntsz = 19;
lnwdth = 1.5;
errWdth = 500; % errorbar width setting

saveFile = ['T:\IDS\Analysis Repository\' num2str(in(1).shot)];

%% Build Figures
[h,ax] = fig_settings(plt,in);
figure(h(1)); % make first figure current

% Pre-Allocation
%RMS = zeros(1+ ~isempty(in(1).doubleplot),...
%    (length(dat(1).impacts))/(1+~isempty(in(n).doubleplot)),length(in),200);

for n = 1:length(in)
%% XXXXXXXXXXXXXXX Main Shot Loop XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    display(['Working on shot: ' num2str(in(n).shot)]);
    saveDat(n).title= in(n).AnalysisTitle;
    saveDat(n).shot=in(n).shot;
    clear data zeroline time
    %addpath('T:\IDS\Data Repository');
    
    %if ~(exist('dat','var') && strcmp(dat(1).title,['Shot ' num2str(in(n).shot)])) % speed up runtime if data is already loaded
        %input=load(['T:\IDS\Data Repository\dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
		%input=load(['C:\Users\Rian\Documents\MATLAB\thosematfilestho\dat' num2str(in(n).shot) '10.mat']); % Local HIT-SI Data (for no BD test)
        clear input dat
        if n==1% Special for noBD test
            input=load(['C:\Users\Rian\Documents\MATLAB\thosematfilestho\dat' num2str(in(n).shot) '10_1.mat']); % Local HIT-SI Data (for no BD test)
        elseif n==2
            input=load(['C:\Users\Rian\Documents\MATLAB\thosematfilestho\dat' num2str(in(n).shot) '.mat']); % Local HIT-SI Data (for no BD test)
        end
        dat=input.dat;
        if flipLoImpact ~=0 % flip the lower array about its centeral impact.
            [dat] = flip_lo(dat,flipLoImpact,in,n);
        end
        dat = trimRange(dat, chan_range, plt.Error,timebound.*(1./in(n).timeScale),[]); % for some reason, this wont save to workspace
        assignin('base','dat',dat);
    %end
    
    Itor(:,n) = dat(1).Itor;
    ItorTime(:,n) = dat(1).ItorTime;
    
    % Plot labels
    switch plt.Type
        case 1
        titles = 'Velocities';
        offset = velSpace;
        units = 'km/s';
        sidebar = [num2str(offset) ' km/s per division'];
        if in(n).error
            errorL = dat(in(n).line).velL;
            errorU = dat(in(n).line).velU;
        end
        case 2
        titles = 'Temperatures';
        offset = tempSpace;
        units = 'eV';
        sidebar = [num2str(offset) ' eV per division'];
        if in(n).error
            errorL = dat(in(n).line).tempL;
            errorU = dat(in(n).line).tempU;
        end
        case 3
        titles = 'Intensities';
        sidebar = ['Arb.'];
        offset = intSpace;
        units = 'Arb.';
        if in(n).error
            errorL = dat(in(n).line).intL;
            errorU = dat(in(n).line).intU;
        end
    end
    if~(in(1).error);errorL=NaN;errorU=NaN;end
    
    
  
    % XXXXXXXXXXXXXXXXXXX Generate Data To Plot XXXXXXXXXXXXXXXXXXXXXX


    % PARAM: [IMPACT, OFFSET, AMP, PHASE, FREQ]


    if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
    %% If we want to doubleplot, but only the raw data

        % Find where each fiber bundle begins and ends.
        doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
        doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

        if plt.Type ==1
            data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,doubleplot(1,:))+in(n).velShift;
            data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
               dat(in(n).line).temp(:,doubleplot(2,:))+in(n).velShift;  
        elseif plt.Type==2   
            data(1:length(dat(1).time),:) = dat(in(n).line).temp(:,doubleplot(1,:));
            data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                dat(in(n).line).temp(:,doubleplot(2,:)); 
        elseif plt.Type==3
             data(1:length(dat(1).time),:) = dat(in(n).line).int(:,doubleplot(1,:));
            data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                dat(in(n).line).int(:,doubleplot(2,:));
            data = in(n).intScale.*data;  
        end


    elseif ~isempty(in(1).fftPlot)
    %% Do the Sine Fit (detect doubleplotting inside)
        dat(in(n).line).vel = averageNans(dat(in(n).line).vel)+in(n).velShift; % remove nans
        dat(in(n).line).temp = averageNans(dat(in(n).line).temp); % remove nans
        %dat(in(n).line).int = averageNans(dat(in(n).line).int).*in(n).intScale; % remove nans
        saveDat(n).rawVel = dat(in(n).line).vel;
        saveDat(n).rawTemp = dat(in(n).line).temp;
        
         
        if ~isempty(in(n).doubleplot)
            % Find where each fiber bundle begins and ends.
            doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
            doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

            % initialize Sine_Fit parameters
            if(n==1);param=zeros(length(doubleplot),10,length(in));end
            param(:,1,n) = dat(1).impacts(doubleplot(1,:));
            param(:,6,n) = dat(1).impacts(doubleplot(2,:));
        else
            doubleplot(1,:) = 1:length(dat(1).impacts);
            if(n==1);param=zeros(length(doubleplot),5,length(in));end
            param(:,1,n) = dat(1).impacts(doubleplot(1,:));
        end
        
         % Param: impacts offset amplitude phase, frequency
         %if ~exist('param','var')
         %       param = zeros(length(doubleplot),10,length(in));
         %end
         
        display('Computing FFT');
        if n==1% initialize
        pRel = zeros(length(doubleplot),1+~isempty(in(n).doubleplot));
        dPar = zeros(length(in),size(doubleplot,2),size(doubleplot,1),4);
        data = zeros((1+~isempty(in(n).doubleplot))*length(dat(1).time),length(doubleplot));
        end
        for i = 1:length(doubleplot) % Loop through impacts ( fits both arrays)
            
            % Fit sine Fn to data, also find phase error
            % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            
            [guess(i,:,n,:),param(i,:,n),saveDat,SigDev(n,i,:),RMS(:,i,n,:),...
                RMS_ideal(:,i,n,:),data(:,i),pRel(i,:),dPar(n,i,:,:)] =...
                sine_fit_module(in, doubleplot,dat,n,i,saveDat,plt);
            
        end
        
        % Find the phase of the injectors
        if n==1
            [injParam] = inj_phase(dat,timebound,in,n);
        end

    else
        %% Just plot Lines
        if plt.Type==1
            data = dat(in(n).line).vel + in(n).velShift;   
        elseif plt.Type==2
            data = dat(in(n).line).temp;
        elseif plt.Type==3
            data = in(n).intScale * dat(in(n).line).int;
        end              
    end

    
    if in(n).doubleplot 
    %% Build timebase
         time(1:length(dat(1).time)) = dat(1).time.*in(1).timeScale + in(n).timeShift;
         time(length(dat(1).time)+1:2*length(dat(1).time)) =  time(1:length(dat(1).time));
    else
        time = dat(1).time.*in(1).timeScale + in(n).timeShift;
    end
    
    % Loop through arrays
    for i = 1: 1+~isempty(in(n).doubleplot)
        %% XXXXXXXXXXXXXXX Profile Calculation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if in(n).fftPlot
        %% If FFT occured, get amp, disp, phase
            dataAvg(:,i) = param(:,2+5*(i-1),n);
            dataStd(:,i) = param(:,3+5*(i-1),n);
            dataPhase(:,i) = param(:,4+5*(i-1),n);
            dataDispl(:,i) = param(:,3+5*(i-1),n).*(1./(param(:,5+5*(i-1),n)*2*pi)) .*1e5;
            
            % Correct the phase measurement to account for 2Pi jumps, has
            % optional sanity check output plots.
            [dataPhase] = correct_phase(dataPhase,plt.SanityPhase,n,i,ax,in,param,data,dat(1).impacts);
            
            % Plot Explainatory Reconstruction
            if n==plt.ExplainReconst % only plot the shot we want
                plotChan=9;
                plotReconst(dat,in,param,n,i,ax,h,plotChan);
            end
            
        else % if no FFT
                % find index corresponding to time bounds
                nTimeLim(1) = find(dat(1).time.*in(n).timeScale >= timebound(1), 1);
                nTimeLim(2) = find(dat(1).time.*in(n).timeScale <= timebound(end), 1, 'last');

                % Calculate all data
                for m = 1:size(data, 2);
                    assignin('base','data',data);
                    selection = data((nTimeLim(1):nTimeLim(2))+(i-1)*size(data,1)/2, m);
                    dataAvg(m,i) = nanmean(selection);
                    dataStd(m,i) = nanstd(selection);
                    dataDispl(m,i) = dataStd(m,i).*(1/14500)./(2*pi) .*1e5;

                end
                
                % Plot raw data
                t2(n) = errorbar(dat(1).impacts(1:size(data,2)), dataAvg(:,i),...
                    dataStd(:,i), 'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        end
    end
    
    
    % XXXXXXXXXXXXXXXXX Plotting Profiles XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % Displacement and FFT are essentially the same for HIT-SI and HIT-SI3
    
    % DISPLACEMENT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if isempty(in(1).fftPlot)
        %% If we arent doing the FFT
        figure(h(4));
        hold on;
        for i = 1: 1+in(n).doubleplot
            plot(dat(1).impacts(1:size(data,2)),dataDispl(:,i),'color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        end
        ylabel('Average Radius of Displacement [cm]');set(gca,'ylim',[0,8]);
        figure(h(2));
    else
        if isempty(in(n).doubleplot); RMS = reshape(RMS,1,[],n,200);end
        [saveDat] = plotDisp(dataDispl,ax,h,i,n,in,saveDat,supress.Disp,...
                    doubleplot,RMS,dat,data,pRel,CutPow,plt.includeTemp,plt.Error,lnwdth);
    end
    
    % FFT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if plt.FFT
        % Plot Normalized Fourrier Power: Injector Frequency/Total
        for i=1:1+in(n).doubleplot
             plot(ax(9),dat(1).impacts(1:size(data,2)),100*pRel(:,i),'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        end
         ylabel(ax(9),'[%]');set(ax(9),'ylim',[0,100]);
         set(ax(9),'xlim',xlim); title(ax(9),'Inj. Mode % of Reconstruction');
        plot(ax(9),xlim,[CutPow,CutPow].*100,'--k');
    end

    
    % Plot Flow Profile if no FFT
    if isempty(in(n).fftPlot) && ~isempty(in(n).doubleplot)
         plotFlow(dataStd,dataAvg,supress.Flow,plt.Error,...
    saveDat,ax,h,n,in,plt.Type,dat,doubleplot,data,[],lnwdth);
       
    elseif ~isempty(in(n).fftPlot)
        
        plot(ax(7),xlim,[0,0],'--k'); set(ax(7),'xlim',xlim);
        ylabel(ax(7),'[km/s]'); set(ax(7),'ylim',[-6,6]);set(ax(7),'xticklabel',[]);
        linkaxes([ax(6) ax(7) ax(8) ],'x');
        if plt.includeTemp
            linkaxes([ax(6) ax(7) ax(8) ax(17) ],'x');
        end
    end
    
    % Check if HIT-SI or HIT-SI3
    if in(n).fftPlot & in(n).doubleplot==1
        %% We're Doing HIT-SI3
        
        % FLOW  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        [saveDat] = plotFlow(param,dataAvg,supress.Flow,plt.Error,...
            saveDat,ax,h,n,in,plt.Type,dat,doubleplot,data,RMS,lnwdth);
        
        % PHASE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        % Initial phase manipulations, Correct_Phase should catch these
        if max(max(dataPhase))>2*pi
            disp(['dataPhase > 2Pi, shifting by: ' num2str(-(max(max(dataPhase))-2*pi))]);
            dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
        end
        % can add 2Pi and maintain relationship with A injector
        if mean(mean(dataPhase)) <=-pi
            dataPhase= dataPhase+2*pi;
            disp('Shifting up by 2pi');
        end
        if mean(mean(dataPhase)) >=pi
            dataPhase= dataPhase-2*pi;
            disp('Shifting down by 2pi');
        end
        
        % Special hardcoded phase changes by 2Pi
        [injParam,dataPhase] = phase_settings(dataPhase,n,in,plt.Type,injParam);
        
        % Flip upper array, look at when toroidal drive happens
        if  driveDirection  == 1
                dataPhase(:,1) = dataPhase(:,1)-pi;
                display('WARNING: POSITIVE FIBER DROPPED -PI, DRIVE DIRECTION PHASE');
        end
        
        % Add in overall phase shift (Usually to flip negative shots phases
        % by 2Pi, so phase is in toroidal current direction).
        dataPhase = dataPhase +in(n).phaseShift;
        
        % Plot Phases
        [saveDat,phaseH] = plotPhase(dat,dataPhase,supress.Phase,plt.Error,...
            ax,h,n,in,data,CutPow,pRel,injParam,xlim,SigDev,doubleplot,lnwdth,saveDat);
        
        % TEMP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if plt.includeTemp
            % Remove Datapoints with too High Error
            if any(supress.Temp(n,:,i))
                dat(in(n).line).temp(:,doubleplot(i,supress.Temp(n,:,i)))=NaN;
            end
            if ~plt.Error % Plot with or without Errorbars
                saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                plot(ax(17),dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            else
                saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                saveDat(n).TempError(:,i) = sqrt(nanmean(dat(in(n).line).tempU(:,doubleplot(i,:)).^2 ))/sqrt(length(dat(1).time));
                saveDat(n).LMTempError(:,i) =nanmean(dat(in(n).line).tempU(:,doubleplot(i,:)));
                errorbar(ax(17),dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))), saveDat(n).TempError(:,i),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            end
            ylabel('[eV]');
            set(ax(6),'xticklabel',[]);
        end
        
    elseif in(n).fftPlot
        %% We're Doing HIT-SI (or only one fiber array)
        display('EXECUTING HIT-SI DATA')
        %RMS = reshape(RMS,1,[],length(in),200);
        % FLOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if any(supress.Flow(n,:))
            param(supress.Flow(n,:),2,n)=NaN;
        end
        if(in(n).velShift~=0);display('VELOCITY SHIFT');end
        if ~plt.Error % Plot with or without errorbars
             saveDat(n).Flow = param(:,2,n);
             t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
        else
            saveDat(n).Flow = param(:,2,n);
            
             error=sqrt( nanmean(dat(in(n).line).velU(:,doubleplot(1,:)).^2) + (squeeze(RMS(1,:,n,100))).^2 )/sqrt(length(dat(1).time));
             saveDat(n).FlowError =error;
             t3(n)=errorbar(ax(7),dat(1).impacts(1:size(data,2)),param(:,2,n),error,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
        end
        set(ax(7),'ylim',[-2,15]);

        % PHASE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        hold on;
        if max(max(dataPhase))>2*pi
            dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
        end
        
        % Eliminate data where the injector frequency spectral power was
        % not over the cutoff level
        for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
        
        % Add in phase shift (usually to flip negative shots by pi)
        dataPhase = dataPhase + in(n).phaseShift;
        
        % Hardcoded Phase Settings
        [injParam,dataPhase] = phase_settings(dataPhase,n,in,plt.Type,injParam);
        
        % Suppress data where the error is too high
        if any(supress.Phase(n,:))
            dataPhase(supress.Phase(n,:),1)=NaN;
        end
        
        if ~plt.Error % Plot errobars
            saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
            phaseH(1,1,n)=plot(ax(8),dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
        else
            saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
            saveDat(n).PhaseError(:,1) = SigDev(n,:,1).*180./pi;
            phaseH(1,1,n)=errorbar(ax(8),dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,squeeze(SigDev(n,:)).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
        end
        ylabel(ax(8),'[deg]');set(ax(8),'ylim',[0,400]);set(ax(8),'xticklabel',[]);
        set(ax(8),'xlim',xlim);
        
        % Plot injector Phase
        if n == 1
            saveDat(n).injPhase = injParam(:,3).*180./pi;
            plot(ax(8),[0 60]',[injParam(1,3) injParam(1,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0    0.4470    0.7410]);
            plot(ax(8),[0 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0.8500    0.3250    0.0980]);
        end
        
        % TEMP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if any(supress.Temp(n,:))
            dat(in(n).line).temp(:,supress.Temp(n,:))=NaN;
        end
        
        if plt.includeTemp
            if ~plt.Error
                saveDat(n).Temp = nanmean(dat(in(n).line).temp);
                plot(ax(17),dat(1).impacts(1:size(data,2)),nanmean(dat(in(n).line).temp),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            else
                saveDat(n).Temp = nanmean(dat(in(n).line).temp);

                saveDat(n).LMTempError = nanmean(dat(in(n).line).tempU);
                display(['TEST: ' num2str(size(dat(in(n).line).tempU)) ' + ' num2str(length(dat(1).time))])
                error=sqrt((nanmean(dat(in(n).line).tempU.^2) + std(dat(in(n).line).temp).^2 )/length(dat(1).time));
                saveDat(n).TempError = error;
                errorbar(ax(17),dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp),error,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            end
            ylabel(ax(17),'[eV]');
            set(ax(6),'xticklabel',[]);
        end
        figure(h(1)) % make other current
    end
    
    
    % Plot Data XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    figure(h(1)); % Make current (unclear why this is necessary)
    [saveDat,t] = plot_data(data,offset,dat,ax,time,errorL,errorU,lnwdth,in,...
        n,saveDat,plt.Tor,timebound,fntsz,plt,sidebar,titles);
end


if ~isempty(in(1).fftPlot)
    %% Add Legend to Phase Plot
    %[leg,icon] = legend(ax8,phaseH(1,:,2),{in(1:1:end).legend});
    % Try this one, may be an R>13 thing
    %[leg,icon] = legend(ax(8),squeeze(phaseH(1,:,end)),{'00000000000000000',...
    %    '00000000000000000'});
    [leg,icon] = legend(ax(8),squeeze(phaseH(1,end,:)),{'00000000000000000',...
        '00000000000000000'});
    sgn = ['-','+'];
    for i = 0: length(in)-1
        try %(~plt.Error || in(1).shot < 829499)
            set(icon(length(in)+1+2*i),'Color',in(i+1).color{1});
            set(icon(length(in)+2+2*i),'Color',in(i+1).color{1});
            set(icon,'LineStyle','-');
        catch
            try
                set(icon(length(in)+1+i).Children.Children,'Color',in(i+1).color{1});
                set(icon(i+(floor(length(icon))/2+1)).Children.Children,'LineStyle','-');
            catch % In case of V<R15
                set( get(get( icon(length(in)+1+i),'Children'),'Children'),'Color',in(i+1).color{1});
                set(get(get(icon(i+(floor(length(icon))/2+1)),'Children'),'Children'),'LineStyle','-');
            end
                
        end
        ItrMean = mean(Itor(find(ItorTime(:,n).*in(1).injTimeScale>timebound(1),1):...
            find(ItorTime(:,n).*in(1).injTimeScale>timebound(2),1),i+1));
        legTxt = [num2str(in(i+1).shot) ': ' sgn((sign(ItrMean)>0)+1) num2str(ItrMean,3)];
        try
            icon(i+1).String=legTxt;
        catch % In case of V<R15
            set(icon(i+1),'String',legTxt)
        end
    end
end
figure(h(1));

%% Figure Properties for Averages (Largely Unused)
if plt.Averages
    figure(h(2)) % make current
    
    title(['Mean ' titles ' and Fluctuations'], 'fontsize', fntsz);
    set(ax(2), 'XLim', [dat(1).impacts(1) - 1, dat(1).impacts(end) + 1]);
    set(gca, 'LineWidth', lnwdth);
    set(gca, 'fontsize', fntsz);
    box on;
    grid on;
    xlabel('Major Radius [cm]');
    ylabel(units);
    
    % Legend
    legendText = cell(1, length(in)); % initialize
    legendHands = zeros(1, length(in)); % initialize
    for n = 1:length(in)
        legendText{n} = in(n).legend;
        legendHands(n) = t2(n);
    end
    legend(ax(2), legendHands, 'Location', 'NorthEastOutside',...
        legendText, 'fontsize', fntsz);
end

%% Plot Currents
ax=plot_currents(plt.compactCurrents,ax,h,dat,in,Itor,ItorTime,fntsz,lnwdth,timebound);

%% Saving
if saving
    cd(['T:\IDS\Analysis Repository\']);
    pause(1.5);% necessary to get figure size correct
    saveas(h(1), ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'fig');
    saveas(h(6), ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
    saveas(h(1), ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'bmp');
    saveas(h(6), ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    if in(1).fftPlot
        saveas(h(7), ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
        saveas(h(7), ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    end
end
