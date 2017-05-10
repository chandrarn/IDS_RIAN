% Convert from .Avi to .mat, for cases where transfering from camera cart
% results in strong compression losses.
% Has a nasty habit of keeping frames before the mark in. 

% Phantom SDK 
addpath(genpath('T:\RChandra\Phantom\PhMatlabSDK'));
LoadPhantomLibraries();
RegisterPhantom(true);

shots = [170502004];

addpath('T:\PhantomMovies\');

for shot = shots
    display(['Shot: ' num2str(shot)]);
    %% Load the frames from the .avi
    v=VideoReader(['T:\PhantomMovies\' num2str(shot) 'f.avi']);
    CineArray = zeros(v.Height,v.Width,v.Duration * v.FrameRate);
    k=1;
    while hasFrame(v)
        CineArray(:,:,k)=double(rgb2gray(readFrame(v)));
        k=k+1;
    end
    % Need to reorient dimensions
    CineArray=shiftdim(CineArray,2);
    %% Load the timebase from the .cine
    [HRES, CH] = PhNewCineFromFile(['T:\PhantomMovies\' num2str(shot) '.cine']);
    EndofFrame = libstruct('tagTIME64');
    pEndofFrame = libpointer('tagTIME64', EndofFrame);
    SelectionCode = PhConConst.GAD_TIME;
    DataSize = uint32(32);
    % Get Trigger Time
    Trigger= libstruct('tagTIME64');
    pTrigger=libpointer('tagTIME64',Trigger);
    PhGetCineInfo(CH, PhFileConst.GCI_TRIGTIMEFR, pTrigger);
    Trigger=pTrigger.Value;   
    Trigger = double(int64((Trigger.seconds*10.^8)))+double(int64((Trigger.fractions*10.^8)/2.^32)); % preform subtraction
    % Exposure
    pExpos = libpointer('uint32Ptr',0);
    PhGetCineInfo(CH, PhFileConst.GCI_EXPOSURENS, pExpos);
    
    TimeVector = zeros((v.Duration * v.FrameRate),1);
    
    [ HRES ] = PhGetCineAuxData( CH, 0, SelectionCode, DataSize, pEndofFrame);
    EndofFrame=pEndofFrame.Value;
    T0=double(int64((EndofFrame.seconds*10.^8)));
    for i = 0:(v.Duration * v.FrameRate)-1
        [ HRES ] = PhGetCineAuxData( CH, i, SelectionCode, DataSize, pEndofFrame);

        EndofFrame=pEndofFrame.Value;
        TimeVector(i+1)= double(int64((EndofFrame.seconds*10.^8))) +double(int64((EndofFrame.fractions*10.^8)/2.^32))-Trigger; %the 2^32 is because Phantom is stupid (check the documentation for TIME64)
    end
     TimeVector = TimeVector-T0;
     Exposure=double((pExpos.Value)/2/1000)*100; % Put Exposure in uS
     TimeVector=TimeVector-Exposure;
     TimeVector=TimeVector/100;
     
     figure; plot(TimeVector);
     title(['Shot: ' num2str(shot)]);
     xlabel('Index');
     ylabel('Time [{\mu}S]');
     grid on;
     
     save(['T:\PhantomMovies\Shot ' num2str(shot) 'f.mat'],'CineArray','TimeVector');
     
%         
    
end
