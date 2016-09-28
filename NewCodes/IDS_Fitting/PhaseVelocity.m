%% Calculate average velocity for a given shot, at a given phase
%  Guess at f with fft, and A with max value
%  fminsearch to get approximate fitting values
%  get the velocity and phase at several different points, average.
%% RIAN CHANDRA, MAY 2014

addpath ('T:\IDS\Data Repository');
shots = [12858710];

TimeBounds= [1,2];
TimeWing = 1/20;

numberSamples=5;


for(i = 1:length(shots))
    
    eval(sprintf('load(''dat%i'');', shot(1))); % get the shot
    
    spacing=linspace(TimeBounds(1),TimeBounds(2),numberSamples);
    for(sample=TimeBounds(1):spacing:TimeBounds(2))
        
        minBound=find(dat.iinjxTime>(sample-.000001 - TimeWing));
        maxBound=find(dat.iinjxTime>(sample-.000001 + TimeWing));
        minBound=minBound(1); maxBound=maxBound(1); %get the element of the first values above the min
        % and max threshold: this should give the bounds sample +- time
        % wing
        f=fft(dat.iinjx(minBound:maxBound));
        f=f(1:  ceil(length(f)/2)); % cut off higher half
        f=find(f>max(real(f))-.01)-1;% should give cycles/durration
        %subtract one because fft starts at 1 cycles/period
        % must have real(f) otherwise everything goes to hell
        f=f/(dat.iinjxTime(maxBound)-dat.iinjxTime(minBound));
        % this should be cycles per mS
        
        A=max(dat.iinjx(minBound:maxBound));

        timespan=(dat.iinjxTime(minBound):.001:dat.iinjxTime(maxBound));
        
        sinfunct= @(x) ((x(1)*sin(2*pi*x(2)*timespan+x(3)))-dat.iinjx(minbound:maxbound)).^2;
        % in theory, this function is the fitted function minus the real
        % data, squared. Essentially, it should give the error squared,
        % which we can attempt to minimize. 
        
        