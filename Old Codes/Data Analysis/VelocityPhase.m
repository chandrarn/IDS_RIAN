% Finding IDS velocity by injector phase.
% Rian Chandra, May 2014

% Function to get the average velocity at a specific injector phase.
% In theory, given that the timestamps for the IDS data acquisition are
% synched with the injectors, all equidistantly time spaced data should be
% at the same phase. We'll see about that. This function should be called
% alongside batch correct, and save velocity and phase seperately, for each
% shot.

%MAY 2014: added 14.5kHz functionality: should be able to save multiple
%phases from the same shot, with cooresponding velocities ( 14.5 has ~9
%velocity points per current period ).

%SMR 2014: added better control for phases between 2pi and zero, added auto
%time bound finding functionality, added fitting error detection, this was
%deleted. 

%% NOTE: Mod(2pi) may be eliminating negative sign for phases. Check this.
function [ Phase,Velocity,Impacts ]= VelocityPhase(dat,TimeBounds,samples,Is14,shift, torPlot)

 close all;
        
    
    
    addpath ('T:\IDS\Data Repository');

    addpath('T:\RChandra\A-A-Ron Code\General Matlab\extrema');
    %TimeBounds= [1,2];


    %Trim bad channels ( copied from Compare_Plots_3 )
    %torPlot = 1; % true for toroidal array (fibers 1 to 36)
    % false for poloidal array (fibers 37 to 72)
    % chan_ranget = [8:24]; % toroidal, mohawk port in midplane
    % chan_ranget = [8:28]; % toroidal, mohawk port perp.
    % chan_ranget = [7:29]; % toroidal, 71 degree port
    % TEMP FOR AARON QUAL FILES
    chan_ranget = [8:24];
    % chan_ranget = [10:27]; % toroidal, axial port
    % chan_ranget = 1:30; % NIMROD mohawk

    timeInMs = 0; % displays time in ms
    shiftTime = 0; % shift time axis for plot a [ms]

    chan_rangep = [46:63]; % poloidal
    chan_rangep = [47:61]; % poloidal
    chan_rangep = [50:58]; % spheromak region
    
    if torPlot
        chan_range = chan_ranget;
    else
        chan_range = chan_rangep;
    end
    % Trim all Data for channel range
    cd('T:\IDS\Display');
    dat = trimRange(dat, chan_range);
    assignin('base','dat',dat);
    dat.time = dat.time + shiftTime;
    %dat.ItorTime = dat.ItorTime + shiftTime;

    if timeInMs
        dat.time = dat.time.*1e3;
    end

    SamplesPerInjectorCycle = 1;
    Impacts = dat.impacts;

    %how to guess frequency: yPos=y.*(y>0); Eliminates times when the
    %graph spkies randomly in the negative region
    %sort imax, disregard first entry and last, take two consecutive
    %values, get x(highervalue)-x(lowervalue) and that should be the
    %period in miliseconds.


    %Guessing phase: if the value is closer to zero than to A and goin
    %g down, its sifted by pi, if goin gup, by 3/2pi, if closer to A
    %and going down, by 2/3pi, going up, by 1/3pi, if negative A: going
    %down, 4/3 pi, going up, 5/3. I think. Math is hard
    %take final result and mod (2pi) to normalize


    % note: if 14.5, in theory, if we start from some time point and go
    % forward by sample multiples of 10, the phases at 1, 11, 21, etc,
    % should all line up, and we can get and average.
    if(Is14) % if 14.5, set end time bound to the velocity time at
             % samples*10 points forward from the initial one
        dummy=find(dat.time>TimeBounds(1)-.001);
        TEMP=dat.time(dummy(1));
        % THE TIMES TEN IS ONLY FOR 14.5, WHERE THERES TEN POINTS PER
        % INJ CYCLE
        temp2 = 0;
        i = 1;
        while i<=samples && temp2<=TimeBounds(2)
            temp2 = dat.time(dummy(1)+10*(i)-1);
            i=i+1
        end% THE MINUS ONE IS IMPORTANT: IE: 100-150 is 51 elements: important for linspace
        TimeBounds(2) = dat.time(dummy(1)+10*(i-2)-1);
        VelminBound = dummy(1)
        VelmaxBound = dummy(1)+10*(i-1) -1;
        samples=i-1 % in case it exits because of time bounts, reset samples: important later
        
        NumVelPoints = samples*10 % How big the velocity array should be
        SamplesPerInjectorCycle=10;
    else
        NumVelPoints=samples;
        VelminBound=find(dat.time>((TimeBounds(1))-.00001));
        VelmaxBound=find(dat.time>((TimeBounds(2))-.00001));
        VelminBound=VelminBound(1)
        VelmaxBound=VelmaxBound(1)
        dat.time(VelminBound)
        dat.time(VelmaxBound)
    end
    colorVec = hsv(samples); % need to assign here: samples can go up


    TimeWing = 1/20; %modefies linspace for x. And like everything else
    Coeffecients=[0];
    Phase=[0];
    Velocity=zeros(size(dat.vel,2),NumVelPoints);
    TEMP=size(Velocity)
    counter=1;

    % get timepoints from the velocity time base ( need this because we
    % only have velocity data for specific timepoints, not as many as
    % for current. Also, these times are in S, not mS. Because.
    % Note: is samples are high enough, it is possible that we will try
    % to read more time points from velPoints than there are dat.time
    % points. In this case, the casting as int should round them to the
    % same number, so we'll just do some twice, which should still be
    % ok, because we average at the end.              
    
    % int64 is important: uint8 cant handle numbers above 255
    VelPoints=int64(linspace(VelminBound,VelmaxBound,NumVelPoints))';

    VelPoints=VelPoints(:)+shift; % shift backwards to set start at phase = 0
    %skip timepoints where dat.vel has NaNs

    i=1;
    while(i<=NumVelPoints)
        i;
        if(find(isnan(dat.vel(VelPoints(i),:))))
            VelPoints=VelPoints+1;
            i=i-1; % This should mean that it checks the same VelPoint again, just incase 
            % the new timepoint ALSO has NaNs
        end
        i=i+1;
    end
    VelPoints

    %% Main loop, get phase angle and velocity values at sample point
    %dat.time(VelPoints).*1e3;%to start at phase=0;
    for(point=VelPoints')%VelPoints(1):VelPoints(end)


            % eval(sprintf('load(''dat%i'');', 12858710)); % get the shot


            minBound=find(dat.iinjxTime>((dat.time(point))-.000001 - TimeWing));

            maxBound=find(dat.iinjxTime>((dat.time(point))-.000001 + TimeWing));
            minBound=minBound(1)
            maxBound=maxBound(1)-1
            % we linspace instead of just using the existing points,
            % min to max, because sometimes y is size 501, and x is
            % only of size 500.
            x=linspace(dat.iinjxTime(minBound),dat.iinjxTime(maxBound),length(minBound:maxBound))';
            %%
            %%
            y=dat.iinjx(minBound:maxBound);

            size(x)
            size(y)
            %guess amplitude
            A=max(dat.iinjx(minBound:maxBound));




%                 %guess frequency
%                 %Get peaks, find distance between them
%                 [xmax,imax]=extrema(y);
%                 imax=sort(imax(1:5));
%                 avdif=0;
%                 N=0;
%                 %take the five highest peaks, sort into chronological
%                 %order, add dist to avdif as long as they arent too close
%                 %(sometimes random spikes can happen, throwing off avdif)
%                 for(i=1:4)
%                     if(imax(i+1)-imax(i)>50)
%                         avdif=avdif+(x(imax(i+1))-x(imax(i))); N=N+1;
%                     end
%                 end
%                 avdif=avdif/N;
%                 f=1/avdif

            %FFT TO FIND FREQUENCY
            NFFT = 2^nextpow2(500); % Next power of 2 from length of y
            Y = fft(dat.iinjx(minBound:maxBound),NFFT)/500;
            %5000=1/time between samples
            f = 5000/2*linspace(0,1,NFFT/2+1);

            %Dont even ask me why this works. I have no idea.
            %Check matlab tutorial for fft
            f=f(find(2*abs(Y(1:NFFT/2+1))==max(2*abs(Y(1:NFFT/2+1)))));




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




         %x=1:.01/pi:2*pi;
         %y=.5*sin(3.*x+pi);
         optimset('TolX',1e-10);
         optimset('MaxIter',1e+10);
         [bestcoeffs,fval,exitflag]=fminsearch(@fun,[A f SinePhase],[],x,y);
         yfit=bestcoeffs(1)*sin(bestcoeffs(2)*2*pi.*x+bestcoeffs(3));
         %Now compare y with yfit
         %figure;
         %plot(x,y,x,yfit);
         %hold on
         %plot([dat.time(point);dat.time(point)],[-A;A],'color','red','marker','*','MarkerSize',18);

         %plot(x,10*sin(2*pi*53.*x+0),'color','red');
         %bestcoeffs
         Coefficients(counter,1)=bestcoeffs(1);
         Coefficients(counter,2)=bestcoeffs(2);
         Coefficients(counter,3)=bestcoeffs(3);


         %Calculate phase at time point,using phase and f, mod 2pi
         % In theory, the phase should be the same for equidistant
         % timepoints. Make sure that this is comparable to other data.
         % We may need to only take time points at the origional points
         % that we actually got data at, not any of the fitted times.
         % those might be the ones which are synched with the injectors
         Phase(counter)=mod(bestcoeffs(2)*2*pi*dat.time(point)+bestcoeffs(3),2*pi);
         Phase(counter);
         %plot([dat.time(point)-TimeWing;dat.time(point)+TimeWing],[bestcoeffs(1)*sin(Phase(counter));bestcoeffs(1)*sin(Phase(counter))],'color','red','marker','o','MarkerSize',18);
        %% CHECK THIS. Should vel for all channels at time point



        Velocity(:,counter)=dat.vel(point,:)';
        %sum(VelocityAvg(14,:));
         %update counter
         counter=counter+1;

    end
    Phase';
    %In theory: this should loop over the 'samples' number of sets of
    %ten phases and velocities, add the 1,11,21,31, etc elemens, and
    %then average them.
    %point
    if(Is14)
        PhaseFinal=zeros(SamplesPerInjectorCycle,samples);%take the median of the  phase values
                                % this takes care of the limiting case
                                % where the phase is ALMOST zero for
                                % one or two data sets, and then 2pi
                                % radians for the other, and the
                                % average is then meaningless.
                                % also may take care of data drift.
        %this takes care of cases where we have less channels
        chan_range(2)-chan_range(1);
        VelocityFinal=zeros(size(dat.vel,2),SamplesPerInjectorCycle);
        size(Velocity);
        size(VelocityFinal);
        Temp=zeros(SamplesPerInjectorCycle,1);
        for(i=0:samples-1)%loop through the 'samples' number of sets
            Temp=Phase((1:SamplesPerInjectorCycle)+SamplesPerInjectorCycle*i)';
            for(j=1:SamplesPerInjectorCycle) % within one set, loop through the ten phases
                Velocity(:,1);

                PhaseFinal(j,i+1)=PhaseFinal(j,i+1)+Phase((i*SamplesPerInjectorCycle)+j);
                (i*SamplesPerInjectorCycle)+j;
                Velocity(:,(i*SamplesPerInjectorCycle)+j);
                size(VelocityFinal(:,j));
                size(Velocity(:,(i*SamplesPerInjectorCycle)+j));
                size(VelocityFinal(:,j));
                VelocityFinal(:,j)=VelocityFinal(:,j)+Velocity(:,(i*SamplesPerInjectorCycle)+j);
            end

        end

        %temporary addition to check for velocity by impact at a given
        %phase
        FVelocity=zeros(size(dat.vel,2),SamplesPerInjectorCycle);
        for(i=1:SamplesPerInjectorCycle)
            TempVel=zeros(size(dat.vel,2),samples);
            size(Velocity);
            for(j=1:samples)
                TempVel(:,j)=Velocity(:,i+SamplesPerInjectorCycle*(j-1))-20;
            end
            STD=std(TempVel');
            MEAN=mean(TempVel,2);
            PHASE=median(PhaseFinal(i,:));
            figure;
            errorbar(Impacts,MEAN',STD);
            xlabel('Impact Parameters','FontSize',15);
            ylabel('Velocity, Km/S, with 1 Std','FontSize',15);
            title1=strcat('Phase: ',num2str(PHASE.*180./(pi)),char(176));
            title(strcat(title1, ' mean velocity'),'FontSize',15);
            figure; hold on;
            size(median(TempVel,2));
            size(FVelocity);
            FVelocity(:,i)=median(TempVel,2);
            for j = 1:samples  
                plot(int64(Impacts),TempVel(:,j),'color',colorVec(j,:));
            end
            plot(int64(Impacts),median(TempVel,2),'-o');
            hold off;
            
            legend('First','Second','Third','Fourth','Fifth','Median');
            title([title1 num2str(samples) ' time points'],'FontSize',15);
            ylabel('Velocity, Km/S','FontSize',15);
            xlabel('Impact Parameters','FontSize',15);

        end

        % There is a particularly bullshit case in which the number of
        % samples is even, and it exactly splits between phases around 6
        % and around zero, so it averages the middle two and gets ~1 pi
        % radians
        if mod(samples,2)==0
            if sum(max(PhaseFinal)>6)==samples/2
                [row,col] = find(PhaseFinal>6,1);
                PhaseFinal(row,col)=0;% This isnt quite science.
            end
        end

        PhaseFinal=median(PhaseFinal,2);%take median of three sets
        Phase=zeros(10);
        %Velocity=zeros(size(dat.vel),samples);
        Phase=(PhaseFinal./(pi));
        %Velocity=VelocityFinal./samples;
        Velocity=FVelocity;
         %order by phase
%         [A,I]=sort(Phase);
%         Phase = Phase(I);
%         Velocity = Velocity(:,I);
    else
        TEMPV=Velocity;
        clear title;
        figure;
        hold on;
        temp=median(Velocity,2);
        plot(linspace(-23,28,size(dat.vel,2)),temp(end:-1:1),'red');
        clear temp;
        temp=mean(Velocity,2);
        plot(linspace(-23,28,size(dat.vel,2)),temp(end:-1:1),'green');
        for i=1:samples
           plot(linspace(-23,28,size(dat.vel,2)),Velocity(end:-1:1,i));
        end
        legend('Median','Mean',strcat('Time Points:  ',num2str(dat.time(VelminBound)),' - ',num2str(dat.time(VelmaxBound))));
        xlabel('Impacts');
        ylabel('Km/S');
        title(dat.title);
        TEMPP=Phase;

        Velocity=median(Velocity,2);%./samples;
        Phase;
        Phase=(sum(Phase)/(samples*pi))'
    end
     %plotting for saftey
     figure;
     plot(Phase);
     title('Phases found')
     ylabel('Multiples of \pi')
     xlabel('Index')
     pause(2); % I want to see the phase graph
     %HOW TO PLOT VELOCITY?
         
        
         
         
  
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



