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

%% XXXXXXXXXXXXXXX Main Shot Loop XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

for n = 1:length(in)
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
    
    %% Plot labels
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
    
    
  
    %% XXXXXXXXXXXXXXXXXXX Generate Data To Plot XXXXXXXXXXXXXXXXXXXXXX


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

    %% Build timebase
    if in(n).doubleplot
         time(1:length(dat(1).time)) = dat(1).time.*in(1).timeScale + in(n).timeShift;
         time(length(dat(1).time)+1:2*length(dat(1).time)) =  time(1:length(dat(1).time));
%              dat(1).time.*in(1).timeScale + in(n).timeShift;
    else
        time = dat(1).time.*in(1).timeScale + in(n).timeShift;
    end

    %% Multiplot has option for plotting data averages by impact





    %% XXXXXXXXXXXXXXX Profile Calculation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    % Loop through arrays
    for i = 1: 1+in(n).doubleplot 
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
            end

end