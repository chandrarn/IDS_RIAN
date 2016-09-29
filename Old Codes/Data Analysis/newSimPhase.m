% Simulation Phase Updated:
% The goal of this progam is to match phases from the simulation data to
% the cooresponding phases in the IDS data for a given shot. It will loop
% across the velocity timebase, find the cooresponding injector time, find
% the cooresponding injector phase a la velocityPhase, do this for the
% whole shot. Then, loop through the phases of the real shot, for each one,
% loop back through the phases that were just found, find the surrounding
% phases, linspace everywhere to get the simulation data for that phase.

% ADD A ZERO TO THE END OF THE SAVE FILE IF CORRECTING 129499 POLOIDAL
function PhaseVelocity = newSimPhase
    torPlot = 0;
    % chan_ranget = 10:35;
    chan_ranget = [8:24]; % toroidal, mohawk port in midplane
    % chan_ranget = [8:28]; % toroidal, mohawk port perp.
    % chan_ranget = [8:27]; % toroidal, 71 degree port
    % chan_ranget = [8:24]; % toroidal, axial port
    % chan_ranget = 1:30; % NIMROD mohawk

    %chan_rangep = 37:72;
    % chan_rangep = [46:63]; % poloidal
    % chan_rangep = [47:58]; % poloidal
     chan_rangep = [50:58]; %poloidal zoomed in on spheromak zone

    close all;
    simDat = 812981010;
    realDat = 129810;
    
    if torPlot
        chan_range = chan_ranget;
    else
        chan_range = chan_rangep;
    end
    
    
    
    % Find simulation Phases
    S = load(['T:\IDS\Data Repository\dat' num2str(simDat) '.mat']);
    dat = S.dat;
    cd('T:\IDS\Display');
    dat = trimRange(dat, chan_range);
    assignin('base','dat',dat)
    TimeWing = .0272;% FOR TET.035; %modefies linspace for x. And like everything else
    for i = 1:length(dat.time)
%         close all;
        %find the indexes of the injector around each timepoint
            i;
            minBound=find(dat.iinjxTime>((dat.time(i))-.000001 - TimeWing));

            maxBound=find(dat.iinjxTime>((dat.time(i))-.000001 + TimeWing));
            minBound=minBound(1);
            maxBound=maxBound(1)-1;

            % This doesnt work in velocityPhase sometimes:
            x = dat.iinjxTime(minBound:maxBound);
            y = dat.iinjx(minBound:maxBound);

            %guess amplitude
            A=max(y);


            %FFT TO FIND FREQUENCY
            NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
            Y = fft(dat.iinjx(minBound:maxBound),NFFT)/length(y);
%            5000=1/time between samples
            f = (1/mean(diff(x)))/2*linspace(0,1,NFFT/2+1);
%            Dont even ask me why this works. I have no idea.
%            Check matlab tutorial for fft
            f=f(find(2*abs(Y(1:NFFT/2+1))==max(2*abs(Y(1:NFFT/2+1)))));
%           Note: the fft program is general, but doesnt always work great.
%           if it fails, just hardcode the frequency in
%              f = 14; % 14kHz


             % Guess phase, taken from VelocityPhase
              %guess phase. In theory, this should produce the correct phase
                %within +- pi/8.  UnitCircle find the (pi/4) section of the unit
                %circle that the first point is in. Going up decides what half of
                %the circle we're on. Phase should correctly add or subract to the
                %correct value, starting at the bottom, 3pi/2.
                %the mod is there because this system does not reset at 2pi, so a
                %y(1) at pi/2 would read 20pi/8

                if(y(1)<-(sqrt(2)/2)*A)
                    UnitCircle=0;
                elseif(y(1)<0)
                    UnitCircle=1;
                elseif(y(1)<(sqrt(2)/2)*A)
                    UnitCircle=2;
                else
                    UnitCircle=3;
                end

                if(y(1)<y(3)) 
                    GoingUp=1;
                else
                    GoingUp=-1;
                end
                
                
                SinePhase = mod((12*pi/8)+ (GoingUp*pi/8) +(GoingUp*UnitCircle*pi/4),2*pi);
                
                InitGuess(i,:) = [A f SinePhase];
                assignin('base','InitGuess',InitGuess);
             
                optimset('TolX',1e-10);
             optimset('MaxIter',1e+10);
             [bestcoeffs,fval,exitflag]=fminsearch(@fun,[A f SinePhase],[],x,y);
             yfit=bestcoeffs(1)*sin(bestcoeffs(2)*2*pi.*x+bestcoeffs(3));
             
%              figure;
%              plot(x,y,'b-o',x,yfit,'r-*',x,(A*sin(f*2*pi.*x + SinePhase)),'g--');
%              hold on
%              plot([dat.time(i);dat.time(i)],[-A;A],'color','red','marker','*','MarkerSize',18);


             % put the cooefecients into an array
             Coefficients(i,1:3)=bestcoeffs(1:3);
             assignin('base','Coefficients',Coefficients);
             Phase(i)=mod(bestcoeffs(2)*2*pi*dat.time(i)+bestcoeffs(3),2*pi);
             Phase(i) = Phase(i)/pi; % important: VelocityPhase does this
             
%              if Phase(i)<.1
%                  figure;
%                  plot(x,y,'b-o',x,yfit,'r-*',x,(A*sin(f*2*pi.*x + SinePhase)),'g--');
%                  hold on
%                  plot([dat.time(i);dat.time(i)],[-A;A],'color','red','marker','*','MarkerSize',18);
%                  Phase(i)
%              end
            
             % %          
            %pause(.5)
    end
    assignin('base','Phase',Phase);
    %figure; plot(Phase);
    
    
    % match sim phases with real phases
    S = load(['T:\IDS\Analysis Repository\Phase Data\Phase' num2str(realDat) '.mat']);
    assignin('base','S',S);
    PhaseVelocity = S.PhaseVelocity;
    finalVelocity = ones(size(dat.vel,2),10,3).*NaN;
    finalSTD = ones(size(dat.vel,2),10,3).*NaN;
    finalPhase = ones(10,3).*NaN;
    for i = 1:length(PhaseVelocity.Phase)
        counter = 1;
        i
        for j = 4:length(dat.time)-5 % cycle through sim phases. 
            if ((Phase(j) < PhaseVelocity.Phase(i)) && (Phase(j+1) > PhaseVelocity.Phase(i))...
                    || ( (Phase(j+1)<Phase(j)) && ( (Phase(j)>PhaseVelocity.Phase(i)) && (Phase(j+1) >PhaseVelocity.Phase(i)))))
                %boolean note: the check is: if we're going up, and the
                %phase we want is inbetween the current phase, and the
                %next one, use the closer of the two. Special case: if you
                %hit the top and go back to the bottom, but your lowest
                %phase is greater than the wanted phase, still trigger.

                % pick the closer of the two phases
                [~,I] = min([(PhaseVelocity.Phase(i)-Phase(j)).^2;(PhaseVelocity.Phase(i)-Phase(j+1)).^2]);
                if I==1
                    finalPhase(i,counter) = Phase(j);% save the phase and velocity just before
                    finalVelocity(:,i,counter) = mean(dat.vel(j-3:j+3,:))';% the IDS phase we want
                    finalSTD(:,i,counter) = std(dat.vel(j-3:j+3,:))';
                    %figure; plot(dat.vel(j-3:j+3,:)'); title(['Phase: ' num2str(Phase(j)) '. TimePoint: ' num2str(j)]);
                else
                    finalPhase(i,counter) = Phase(j+1);% save the phase and velocity just before
                    finalVelocity(:,i,counter) = mean(dat.vel(j-2:j+4,:))';% the IDS phase we want
                    finalSTD(:,i,counter) = std(dat.vel(j-2:j+4,:))';
                    %figure; plot(dat.vel(j-3:j+3,:)'); title(['Phase: ' num2str(Phase(j)) '. TimePoint: ' num2str(j)]);
                end
                 counter = counter + 1;
%                  if i ==7 % temp test to find the y-inj
%                      figure; plot(dat.vel(j-3:j+3,:));
%                      figure; plot(dat.vel(j-10:j+10,:));
%                  end
%                  if j >40 && j<60
%                      figure; plot(dat.vel(j-10:j+10,:));
%                      title(num2str(Phase(j)));
%                  end
            end
        end
        figure; errorbar([squeeze(finalVelocity(:,i,:)),median(finalVelocity(:,i,:),3)],[squeeze(finalSTD(:,i,:)),std(finalVelocity(:,i,:),[],3)]);
        title(['Phase: ' num2str(PhaseVelocity.Phase(i))])
        %set(gca,'XLim',[50 58]);
    end
    assignin('base','finalPhase',finalPhase);
    assignin('base','finalVelocity',finalVelocity);
    assignin('base','finalSTD',finalSTD);
    
    % resave the data in a struct
    PhaseVelocity.Phase = median(finalPhase,2);
    PhaseVelocity.Std = std(finalVelocity,[],3);
    PhaseVelocity.Velocity = median(finalVelocity,3);
    
    save(['T:\IDS\Analysis Repository\Phase Data\Phase' num2str(simDat) '.mat'],'PhaseVelocity');
    assignin('base','PhaseVelocity',PhaseVelocity);
    
     % Copied from velocityPhase           
     % actual sine function. Has to be sepearate from fminsearch because stupid
     function out=fun(coeff,X,Y)
         a = coeff(1);
         b = coeff(2);
         c = coeff(3);
         Y_fun = a .* sin(b*2*pi.*X+c);
         DIFF = Y-Y_fun; 
         SQ_DIFF = DIFF.^2;

         out = sum(SQ_DIFF);
     end
    
 
end