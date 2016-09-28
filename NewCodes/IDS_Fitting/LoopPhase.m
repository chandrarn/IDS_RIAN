% loop through shots, get phase and velocity

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
shots = [129810]
samples = 5;
shift = -3;% shift the start time, to start at phase = 0;
Is14=1;
addpath('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes');

    for( i = shots)
        %cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data')
        cd('T:\IDS\Data Repository\TEMP');
        eval(sprintf('load(''dat%i10'');', i));
        [Phase,Velocity,Impacts]=VelocityPhase(dat,[1.49 1.87],samples,Is14,shift);
        cd('T:\RChandra\A-A-Ron Code\Matlab Code\Analysis Codes\Phase Data\Temp');
        PhaseVelocity.Phase=Phase;
        PhaseVelocity.Velocity=Velocity;
        PhaseVelocity.Impacts=Impacts;
        save(['Phase' num2str(i) ], 'PhaseVelocity');
    end


        
        