% Rian Chandra
% loop through various shots, finds the median velocity at a specific
% phase. Saves the phase, the velocity, and the impacts for the shots. Can
% no longer save error. 

% new hitsi3 shots
% 140701033: He, -15 Phase: NO DATA
% 140701034: He, -18 Phase: NO DATA
% 140702015: He, -22 Phase: 0
% 140702026: He, 27  Phase: NO DATA
% 140708016: He, 26  Phase: 0
% 140708017: He, 27  Phase: 0
% 140708018: D   -29 Phase: 0
% 140709010: He  12  Phase: 120
% 140709012: He  20  Phase: 120
% 140709013: He  12  Phase: 120


%shots = [140701033,140701034,140702015,140702026,140708016,140708017,140708018,140709010,140709012,140709013];
%shots = [129495, 129496, 129793,   129440, 129441, 129443, 129446, 129449, 129450, 129451];
%time = [1.57,1.94; 1.51,1.94; 1.4,1.85;   1.74,1.95; 1.54,1.94; 1.34,1.94; 1.61,1.93; 1.66,1.95; 1.68,1.94; 1.41,1.95 ];
%shots = [129518, 129523, 129528, 129530, 129592];
%time = [1.21,1.94; 1.42,1.93; 1.51,1.91; 1.51,1.93; 1.45,1.67;] ;
%shots = [129810,129817, 129819, 129820];
%time = [1.47,1.85; 1.2,1.85; 1.39,1.84; 1.2,1.75 ];
shots = [129530];
time = [1.5,1.9];
samples = 5;
shift = 0;% shift the start time, to start at phase = 0;
Is14=1;
addpath('T:\IDS\Data Analysis\');
torPlot = 0;

    for( i = 1:length(shots))
        display(['Working on shot: ' num2str(shots(i))]);
        %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')
        cd('T:\IDS\Data Repository\');%\TEMP
        eval(sprintf('load(''dat%i10'');', shots(i)));
        [Phase,Velocity,Impacts]=VelocityPhase(dat,time(i,:),samples,Is14,shift,torPlot);
        cd('T:\IDS\Analysis Repository\Phase Data');
        PhaseVelocity.Phase=Phase;
        PhaseVelocity.Velocity=Velocity;
        PhaseVelocity.Impacts=Impacts;
        %PhaseVelocity.Std = Std; % I overwrote std functionality in
        %velocity phase like an idiot
        assignin('base','PhaseVelocity',PhaseVelocity);
        %save(['Phase1' num2str(shots(i)) ], 'PhaseVelocity');
    end


        
        