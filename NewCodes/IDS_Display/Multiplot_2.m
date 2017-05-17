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

%% Main Loop

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
    
    
    if plotType==1 || plotType ==2 || plotType==3
    %% Define the plotting Data
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
                        sine_fit(in, doubleplot,dat,n,i,saveDat);
                
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
            
    end
end