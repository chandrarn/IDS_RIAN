%{ 
MultiPlots Improved:
Modularized the input settings, and fiber flip

Directly connected to IAEA Plots


%}

%% Settings
% Pick set of shots and analysis settings
 n = 1;% 160728013, 160728012
% n = 7;% 129499, 129496
[in,timebound,chan_range,xlim,supress] = in_settings(n);

% Pick plot type
plt.Type = 1; % Velocity
% plotType = 2; % Temperature
% plotType = 3; % Intensity
%plotType = 4; % Displacement

% Which plots to make, other settings
saving = 0;
plt.FFT = 0;
plt.Currents = 1;
plt.SanityPhase = 0;
plt.Averages = 0;
plt.compactCurrents = 1;
plt.Tor = 0; % The line plot will show the differenve between fibers
plt.ExplainReconst=0;% will plot the data and reconstruction from line n to a seperate plot
driveDirection=1*(pltType==1); % upper fiber dropped by -pi, phase now references "positive toroidal flow magnitude"
flipLoImpact=[47].*(in(1).shot>229499); % correct for lower fiber being rotated 180 about impact 30 
% the bracketed number is the pre-trimmed index of this impact. The
% following data is hardcoded for 160728017 data range
plt.Error = 1;
plt.includeTemp = 1; % include temperature in the velocity plot
CutPow  = .1; % FFT Power Cuttoff Percentage. 

if isempty(in(1).fftPlot)
    Analysis=1;
else
    Analysis = 2; % Analyze torroidal flow, Amplitude, phasing. Replaces "Averages"
end

% Spacing for line plot
velSpace = 15; % km/s
intSpace = 20; % arb.u.
tempSpace = 20; % eV

saveFile = ['T:\IDS\Analysis Repository\' num2str(in(1).shot)];

%% Build Figures
[h,ax] = fig_settings(plt,Analysis,in)
figure(h(1)); % make first figure current



for n = 1:length(in)
%% XXXXXXXXXXXXXXX Main Shot Loop XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    n
    saveDat(n).title= in(n).AnalysisTitle;
    saveDat(n).shot=in(n).shot;
    clear data zeroline time
    %addpath('T:\IDS\Data Repository');
    
    if ~(exist('dat','var') && strcmp(dat(1).title,['Shot ' num2str(in(n).shot)])) % speed up runtime if data is already loaded
        input=load(['T:\IDS\Data Repository\dat' num2str(in(n).shot) '10.mat']); % Real HIT-SI Data
        dat=input.dat;
        if flipLoImpact ~=0 % flip the lower array about its centeral impact.
            [dat] = flipLo(dat,flipLoImpact);
        end
        dat = trimRange(dat, chan_range, plotError,timebound.*(1./in(n).timeScale),[]); % for some reason, this wont save to workspace
        assignin('base','dat',dat);
    end
    
    Itor(:,n) = dat(1).Itor;
    
    % Plot labels
    if plotType==1
        titles = 'Velocities';
        offset = velSpace;
        units = 'km/s';
        sidebar = [num2str(offset) ' km/s per division'];
        if in(n).error
            errorL = dat(in(n).line).velL;
            errorU = dat(in(n).line).velU;
        end
    elseif plotType==2
        titles = 'Temperatures';
        offset = tempSpace;
        units = 'eV';
        sidebar = [num2str(offset) ' eV per division'];
        if in(n).error
            errorL = dat(in(n).line).tempL;
            errorU = dat(in(n).line).tempU;
        end
    elseif plotType==3
        titles = 'Intensities';
        sidebar = ['Arb.'];
        offset = intSpace;
        units = 'Arb.';
        if in(n).error
            errorL = dat(in(n).line).intL;
            errorU = dat(in(n).line).intU;
        end
    end
    
    
  
    % XXXXXXXXXXXXXXXXXXX Generate Data To Plot XXXXXXXXXXXXXXXXXXXXXX


    % PARAM: [IMPACT, OFFSET, AMP, PHASE, FREQ]


    if ~isempty(in(n).doubleplot) && isempty(in(1).fftPlot)
    %% If we want to doubleplot, but only the raw data

        %data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,in(n).doubleplot(1,:));
        %data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
        %    dat(in(n).line).vel(:,in(n).doubleplot(2,:));
        doubleplot(1,:) = 1:(length(dat(1).impacts))/2;
        doubleplot(2,:) = (length(dat(1).impacts)/2)+1:length(dat(1).impacts);

        if plotType ==1
            data(1:length(dat(1).time),:) = dat(in(n).line).vel(:,doubleplot(1,:))+in(n).velShift;
            data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
               dat(in(n).line).temp(:,doubleplot(2,:))+in(n).velShift;  
        elseif plotType==2   
            data(1:length(dat(1).time),:) = dat(in(n).line).temp(:,doubleplot(1,:));
            data(length(dat(1).time)+1:2*length(dat(1).time),:) = ...
                dat(in(n).line).temp(:,doubleplot(2,:)); 
        elseif plotType==3
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
            param(:,1,n) = dat(1).impacts(doubleplot(1,:));
            param(:,6,n) = dat(1).impacts(doubleplot(2,:));
        else
            doubleplot(1,:) = 1:length(dat(1).impacts);
            param(:,1,n) = dat(1).impacts(doubleplot(1,:));
        end
%                 if flipLoImpact == 1 % flip the lower array about its centeral impact.
%                     dat(in(n).line).vel(:,doubleplot(2,:)) = dat(in(n).line).vel(:,doubleplot(2,end):-1:doubleplot(2,1));
%                     dat(in(n).line).temp(:,doubleplot(2,:)) = dat(in(n).line).temp(:,doubleplot(2,end):-1:doubleplot(2,1));
%                 end
         % Param: impacts offset amplitude phase, frequency
         if ~exist('param','var')
                param = zeros(length(doubleplot),10,length(in));
         end
        display('Computing FFT');
        pRel = zeros(length(doubleplot),2);
        dPar = zeros(length(in),size(doubleplot,2),size(doubleplot,2),4);

        for i = 1:length(doubleplot) % Loop through impacts ( fits both arrays)

            [guess(i,:,n,:),param(i,:,n),saveDat,SigDev(n,i,:),RMS(:,i,n,:),RMS_ideal(:,i,n,:)] =...
                sine_fit_module(in, doubleplot,dat,n,i,saveDat);

        end

        [injParam] = inj_phase(dat,timebound,in);

    else
        %% Just plot Lines
        if plotType==1
            data = dat(in(n).line).vel + in(n).velShift;   
        elseif plotType==2
            data = dat(in(n).line).temp;
        elseif plotType==3
            data = in(n).intScale * dat(in(n).line).int;
        end              
    end

    
    if in(n).doubleplot 
    %% Build timebase
         time(1:length(dat(1).time)) = dat(1).time.*in(1).timeScale + in(n).timeShift;
         time(length(dat(1).time)+1:2*length(dat(1).time)) =  time(1:length(dat(1).time));
%              dat(1).time.*in(1).timeScale + in(n).timeShift;
    else
        time = dat(1).time.*in(1).timeScale + in(n).timeShift;
    end

    
    % Multiplot has option for plotting data averages by impact

    
    % Loop through arrays
    for i = 1: 1+in(n).doubleplot
        %% XXXXXXXXXXXXXXX Profile Calculation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if in(n).fftPlot
        %% If FFT occured, get amp, disp, phase
            dataAvg(:,i) = param(:,2+5*(i-1),n);
            dataStd(:,i) = param(:,3+5*(i-1),n);
            dataPhase(:,i) = param(:,4+5*(i-1),n);
            dataDispl(:,i) = param(:,3+5*(i-1),n).*(1./(param(:,5+5*(i-1),n)*2*pi)) .*1e5;
            
            % Correct the phase measurement to account for 2Pi jumps, has
            % optional sanity check output plots.
            [dataPhase] = correct_phase(dataPhase,plotSanityPhase,n,i,ax,in,param);
            
            % Plot Explainatory Reconstruction
            if n==plotExplainReconst % only plot the shot we want
                plotChan=9;
                plotReconst(dat,in,param,n,i,ax,h,plotChan);
            end
            
        else % if no FFT
                % find index corresponding to time bounds
                nTimeLim(1) = find(dat(1).time.*in(n).timeScale >= timebound(1), 1);
                nTimeLim(2) = find(dat(1).time.*in(n).timeScale <= timebound(end), 1, 'last');
                % nTimeLim(2) = nTimeLim(2) - 1; % the above command goes one too far

                % Calculate all data
                for m = 1:size(data, 2);
                    assignin('base','data',data);
                    selection = data((nTimeLim(1):nTimeLim(2))+(i-1)*size(data,1)/2, m);
                    dataAvg(m,i) = mean(selection(~isnan(selection)));
                    dataStd(m,i) = std(selection(~isnan(selection)));
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
    [saveDat] = plotDisp(dataDispl,ax,h,i,n,in,saveDat,dispSupress,...
                doubleplot,RMS,dat,data,pRel,CutPow,includeTemp);
    
    % FFT XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if plotFFT
        for i=1:1+in(n).doubleplot
             plot(ax(9),dat(1).impacts(1:size(data,2)),100*pRel(:,i),'-*','color', in(n).color{i}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
        end
         ylabel(ax(9),'[%]');set(ax(9),'ylim',[0,100]);
         set(ax(9),'xlim',xlim); title(ax(9),'Inj. Mode % of Reconstruction');
        plot(ax(9),xlim,[CutPow,CutPow].*100,'--k');
    end
    if n==1;figure(h2);end
    
    % Plot Flow Profile if no FFT
    if isempty(in(n).fftPlot) && in(n).doubleplot==1

        plot(dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2))./2,'color' ,'k','marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});
        plot(xlim,[0,0],'--k')
        ylabel('Toroidal Flow [km/s]'); set(gca,'ylim',[-10,10]);
        xlabel('Impacts [cm]');
        
    % Plot the black dashed line on the flow plot otherwise
    elseif ~isempty(in(n).fftPlot)
        
       % plot(ax7,dat(1).impacts(1:size(data,2)),-(dataAvg(:,1)-dataAvg(:,2)),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style);
        plot(ax7,xlim,[0,0],'--k'); set(ax7,'xlim',xlim);
        ylabel(ax7,'[km/s]'); set(ax7,'ylim',[-6,6]);set(ax7,'xticklabel',[]);
        h2.delete;
        linkaxes([ax6 ax7 ax8 ],'x');
        if includeTemp
            linkaxes([ax6 ax7 ax8 ax17 ],'x');
        end
    end
    
    % Check if HIT-SI or HIT-SI3
    if in(n).fftPlot & in(n).doubleplot==1
        %% We're Doing HIT-SI3
        
        % FLOW  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        [saveDat] = plotDisp(dataDispl,ax,h,n,in,saveDat,dispSupress,...
            doubleplot,RMS,dat,data,pRel,CutPow,includeTemp);
        
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
        [injParam,dataPhase] = phase_settings(dataPhase,n,in,plotType,injParam);
        
        % Flip upper array, look at when toroidal drive happens
        if  driveDirection  == 1
                dataPhase(:,1) = dataPhase(:,1)-pi;
                display('WARNING: POSITIVE FIBER DROPPED -PI, DRIVE DIRECTION PHASE');
        end
        dataPhase = dataPhase +in(n).phaseShift;
        
        % Plot Phases
        [saveDat,phaseH] = plotPhase(dat,dataPhase,phaseSupress,plotError,...
            ax,h,n,data,CutPow,pRel,injParam,xlim,SigDev);
        
        % TEMP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if includeTemp
            if any(tempSupress(n,:,i))
                dat(in(n).line).temp(:,doubleplot(i,tempSupress(n,:,i)))=NaN;
            end
            if ~plotError
                saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                plot(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            else
                saveDat(n).Temp(:,i) = mean(dat(in(n).line).temp(:,doubleplot(i,:)));
                saveDat(n).TempError(:,i) = sqrt(nanmean(dat(in(n).line).tempU(:,doubleplot(i,:)).^2 ))/sqrt(length(dat(1).time));
                saveDat(n).LMTempError(:,i) =nanmean(dat(in(n).line).tempU(:,doubleplot(i,:)));
                errorbar(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp(:,doubleplot(i,:))), saveDat(n).TempError(:,i),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            end
            ylabel('[eV]');
            set(ax6,'xticklabel',[]);
        end
        
    elseif in(n).fftPlot
        %% We're Doing HIT-SI
        display('EXECUTING HITSI DATA')
        
        % FLOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if any(flowSupress(n,:))
            param(flowSupress(n,:),2,n)=NaN;
        end
        if ~plotError
             saveDat(n).Flow = param(:,2,n);
             t3(n)=plot(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
        else
            saveDat(n).Flow = param(:,2,n);
             error=sqrt( nanmean(dat(in(n).line).velU(:,doubleplot(1,:)).^2) + (squeeze(RMS(1,i,n,100))).^2 )/sqrt(length(dat(1).time));
             saveDat(n).FlowError =error;
             t3(n)=errorbar(ax7,dat(1).impacts(1:size(data,2)),param(:,2,n),error,'color',[in(n).color{1}],'marker','*','LineWidth', lnwdth, 'LineStyle', in(n).style{1});%,...
        end
        set(ax7,'ylim',[-2,15]);

        % PHASE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        hold on;
        if max(max(dataPhase))>2*pi
            dataPhase=dataPhase-(max(max(dataPhase))-2*pi);
        end
        for j=1:length(dataPhase(:,1));try dataPhase(j.*(pRel(j,i)<CutPow),1)=NaN;end;end
        dataPhase = dataPhase + in(n).phaseShift;
        % Hardcoded Phase Settings
        [injParam,dataPhase] = phase_settings(dataPhase,n,in,plotType,injParam);
        
        if any(phaseSupress(n,:))
            dataPhase(phaseSupress(n,:),1)=NaN;
        end
        if ~plotError
            saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
            phaseH(1,1,n)=plot(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
        else
            saveDat(n).Phase(:,1) = dataPhase(:,1).*180./pi;
            saveDat(n).PhaseError(:,1) = SigDev(n,:,1).*180./pi;
            phaseH(1,1,n)=errorbar(ax8,dat(1).impacts(1:size(data,2)),dataPhase(:,1).*180./pi,squeeze(SigDev(n,:)).*180./pi,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{1});
        end
        ylabel(ax8,'[deg]');set(ax8,'ylim',[0,400]);set(ax8,'xticklabel',[]);
        set(ax8,'xlim',xlim);
        % Plot injector Phase
        if n == 1
            saveDat(n).injPhase = injParam(:,3).*180./pi;
            plot(ax8,[0 60]',[injParam(1,3) injParam(1,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0    0.4470    0.7410]);
            plot(ax8,[0 60]',[injParam(2,3) injParam(2,3)]'.*180./pi,'--','LineWidth',lnwdth,'color',[    0.8500    0.3250    0.0980]);
        end
        
        % TEMP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        if any(tempSupress(n,:))
            dat(in(n).line).temp(:,tempSupress(n,:))=NaN;
        end
        if includeTemp
            if ~plotError
                saveDat(n).Temp = nanmean(dat(in(n).line).temp);
                plot(ax17,dat(1).impacts(1:size(data,2)),nanmean(dat(in(n).line).temp),'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            else
                saveDat(n).Temp = nanmean(dat(in(n).line).temp);

                saveDat(n).LMTempError = nanmean(dat(in(n).line).tempU);
                display(['TEST: ' num2str(size(dat(in(n).line).tempU)) ' + ' num2str(length(dat(1).time))])
                error=sqrt((nanmean(dat(in(n).line).tempU.^2) + std(dat(in(n).line).temp).^2 )/length(dat(1).time));
                saveDat(n).TempError = error;
                %sqrt(mean(dat(in(n).line).tempU(:,doubleplot(i,:)).^2 ))/sqrt(length(dat(1).time));
                errorbar(ax17,dat(1).impacts(1:size(data,2)),mean(dat(in(n).line).temp),error,'-*','color', in(n).color{1}, 'LineWidth', lnwdth, 'LineStyle', in(n).style{i});
            end
            ylabel(ax17,'[eV]');
            set(ax6,'xticklabel',[]);
        end
        figure(h) % make other current
    end
    
    % Plot Data
    [saveDat,t] = plot_data(data,offset,dat,ax,time,errorL,errorU,lnwdth,in,n,saveDat,t,plotTor);
end


if ~isempty(in(1).fftPlot)
    %% Add Legend to Phase Plot
    %[leg,icon] = legend(ax8,phaseH(1,:,2),{in(1:1:end).legend});
    [leg,icon] = legend(ax8,phaseH(1,1,:),{in(1:1:end).legend});
    for i = 0: length(in)-1
        try %(~plotError || in(1).shot < 829499)
            set(icon(length(in)+1+2*i),'Color',in(i+1).color{1});
            set(icon(length(in)+2+2*i),'Color',in(i+1).color{1});
            set(icon,'LineStyle','-');
        catch
            set(icon(length(in)+1+i).Children.Children,'Color',in(i+1).color{1});
            set(icon(i+(floor(length(icon))/2+1)).Children.Children,'LineStyle','-');
        end
    end
end
figure(h(1));

%% Figure Properties for Averages
if plotAverages
    figure(h2) % make current
    
    title(['Mean ' titles ' and Fluctuations'], 'fontsize', fntsz);
    set(ax2, 'XLim', [dat(1).impacts(1) - 1, dat(1).impacts(end) + 1]);
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
    legend(ax2, legendHands, 'Location', 'NorthEastOutside',...
        legendText, 'fontsize', fntsz);
end


%% Saving
cd(['T:\IDS\Analysis Repository\']);
if saving
    pause(1.5);% necessary to get figure size correct
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'fig');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
    saveas(h, ['Lines L' num2str(in(1).line) titles num2str(in(1).shot)], 'bmp');
    saveas(h6, ['Analysis L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    if in(1).fftPlot
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'fig');
        saveas(h7, ['Spectrum L' num2str(in(1).line) num2str(in(1).shot)], 'bmp');
    end
end
