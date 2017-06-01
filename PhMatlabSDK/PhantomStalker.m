%%%CINE-TO-MAT CONVERTER PROGRAM
%%RIAN CHANDRA 2014 & KYLE ROBERTS 2016 & RIAN CHANDRA 2016
%{
This code converts and saves movies from .cine to .mat file format.
It can handle old files, or pull them from the camera as they are taken.
Sequence of events:
- (if new shots) Setup camera, get shot number
- Make CineArray and TimeVector, reset image processing parameters
- Loop through frames, take the pixel values and time.
- Save to timestamped folder
%}

%% The idea is to have some kind of while loop that allows the program to run in the background and wait for a shot.
clear all;
%close all;
%clc;
import MDSplus.*
addpath(genpath('C:\Users\HITSI\Documents\GitHub\IDS_RIAN\PhMatlabSDK'));
pdc3 =0; % pdc3 flag should be set automatically
% Correcting old data flag
OldData = 0;
shots =  [170519001:170519002];% Old data to correct, in ascending order.
SavePath = 'T:\\PhantomMovies\\';
% Load phaotm librarie
LoadPhantomLibraries();
RegisterPhantom(true);

%% Getting time signatures associated with begining and ending of movie
if ~OldData % Pull the shot length from the tree
    Conn = Connection('landau.hit');
    Conn.openTree('hitsi3',0);
    MovieStop = double(Conn.get('\Shot_Length')) + 0.0005;
    MovieStop = MovieStop*10^(6); % unit conversion
    % Refresh CCD before shot is taken
    HRES=calllib('phcon', 'PhBlackReferenceCI', 0,[]); % CSR
    if HRES == 132
        disp('Successfully refreshed CCD');
    elseif HRES == -259 % This may occur if matlab is left on overnight
        error('No Responce From Camera, Restart Matlab');
    end
    PhRecordCine(0); % RNC said to add this 8/2/16 - will pretrigger camera
    
    % Check PDC3 and HITSI3 shot number ( Wont work before baseline shot is
    % taken)
    try
        Conn.openTree('pdc3',0);
        pdc3Shot = double(Conn.get('$Shot'));
    catch
        pdc3Shot = NaN; % if no PDCbaseline shot, this wont work.
    end
    Conn.openTree('hitsi3',0);
    hitsiShot = double(Conn.get('$Shot'));
        
else
    % If correcting old data, assume that it wont be longer than a minute
    MovieStop = 60*10^6;
end

%% Make the folder to save in 
 % This creates a new folder for the day, for the cines to be saved in
t = datestr(now, 'mmddyy'); %the date, with the format in quotes 
files = dir('E:\\Cines\\'); 
W = ['HS Camera ' t '\\']; %name of the new folder 'HS Camera mmddyy'

t = 0; %create the flag
for i = 1:length(files) %checks to see if that folder already exists
    if strcmp(files(i).name,W)
        t=t+1;
    end
end

if t < 1 %if the folder doesn't exist
    mkdir('E:\\Cines\\',W); %makes a new folder in the cines drive
end
SaveCinePath = ['E:\\Cines\\' W];

%% Main Loop Over Shots, Old or Current
cont = 1;shotnum=0;
emergencycounter = 0;

while cont==1;
    display('Waiting for store to conclude');
    % Wait for the shot to happen
    if ~OldData
        % The tree will send out this flag after a shot in hitsi3 and pdc3
        Conn.openTree('hitsi3',0); 
        Conn.get('wfevent("DTACQ_TRIGGERED")');
        
    end
    
    if OldData
            % Get the next shot and start to convert it
            shotnum=shots(find(shots>=(shotnum+1),1,'first'));
            if shotnum==shots(end) % if we've reached the end, stop
                cont=0;
            end
    else
        % Try to get the new shot number
        try
            Conn.openTree('pdc3',0); % Note: There is a discrepancy between 
            pdcNew = double(Conn.get('$Shot'));% the pdc3 shotnumber and the 
        catch                                  % Python script number
            pdcNew = NaN;
        end
        Conn.openTree('hitsi3',0);
        hitsiNew = double(Conn.get('$Shot'));
        if hitsiNew ~= hitsiShot % check to see if the hitsi3 shot changed
            shotnum = hitsiNew;
            pdc3=0; % set flag for storing to tree (PDC doesnt store IDS to tree)
        elseif ~(isnan(pdcNew)) && (pdcNew ~= pdc3Shot)  % Othewise its pdc3
            shotnum = pdcNew;
            pdc3=1;
        else
            disp('Failed to get new shot number, trying again');
            shotnum = hitsiShot + emergencycounter;
            emergencycounter=emergencycounter+1; % this can be used if you REALLY
            % need to get data and it isn't getting the shot number. 
            %continue; % Skip rest of program, wait for next shot
        end
        hitsiShot=hitsiNew; % update the shot number
        pdc3Shot = pdcNew;
    end
    
    %% Shot happened, begin conversion
    %clearvars -except 'Exposure' 'Shots' 'j' 'ShotNameVector' 'date' ;
    % Figure out how to get the Shot number
    fprintf('Working on Cine #%d...\n', shotnum)
    
    %This creates a new handle for the cine 1, in camera 0. The handle is CH
    if ~OldData % get pointer to new cine in camera
        [HRES, CH] = PhNewCineFromCamera(0,1);
    else % Get pointer to old cine to convert
        [HRES, CH] = PhNewCineFromFile(['T:\PhantomMovies\' num2str(shotnum) '.cine']);
        if isempty(CH)
            display(['Shot: ' num2str(shotnum) ' Not Found']);
            break
        end
    end
    
    % Reset Image Processing Parameters. THIS IS IMPORTANT: if these are
    % not reset, image processing done in the PCC will propagate through to
    % the saved .mat file, rendering it unusable. 
    PhSetCineInfo(CH,PhFileConst.GCI_BRIGHT,libpointer('float',0.0));
    PhSetCineInfo(CH,PhFileConst.GCI_CONTRAST,libpointer('float',1.0));
    PhSetCineInfo(CH,PhFileConst.GCI_GAMMA,libpointer('float',1.0));
    
    %% Setup CineArray    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Converting Frames');
    %Get first frame number
    pFirstIm = libpointer('int32Ptr',0);
    PhGetCineInfo(CH, PhFileConst.GCI_FIRSTIMAGENO, pFirstIm);
    firstIm = pFirstIm.Value;
    %get last frame numberS
    pImCount = libpointer('uint32Ptr',0);
    PhGetCineInfo(CH, PhFileConst.GCI_IMAGECOUNT, pImCount);
    lastIm = int32(double(firstIm) + double(pImCount.Value) - 1);
    %Get frame size
    pHeight = libpointer('uint32Ptr',0);
    pWidth = libpointer('uint32Ptr',0);
    PhGetCineInfo(CH, PhFileConst.GCI_IMHEIGHT,pHeight);
    PhGetCineInfo(CH, PhFileConst.GCI_IMWIDTH, pWidth);
    Height=pHeight.Value;
    Width=pWidth.Value;
    %Create array "CineArray" to put frames into of correct size
    %CineArray =zeros(Height,Width,pImCount.Value);
    CineArray =zeros(Height,Width,lastIm);
    %make the array out of the frames ( origIm gets thrown out)
    origIm=zeros(Height,Width);
    %Get the frame bit info
    sFrame = libstruct('tagIMRANGE');
    sFrame.cnt=1;
    sFrame.First = 0;
    [HRES, pixels, IH] = PhGetCineImage(CH,sFrame,(50000000));
    indNum = IH.biSizeImage; % total number of data blocks in image (2x number of pixels)
    
    
    %% Make time vector    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Creating Timebase');
    %Create variables to select and store the end-of-frame times
    EndofFrame = libstruct('tagTIME64');
    pEndofFrame = libpointer('tagTIME64', EndofFrame);
    SelectionCode = PhConConst.GAD_TIME;
%     ExposureTime = PhConConst.GAD_EXPOSURE;
    DataSize = uint32(32);
    exposure = uint16(16);
    %Get total number of frames
    pImCount = libpointer('uint32Ptr',0);
    pExpos = libpointer('uint32Ptr',0);
    PhGetCineInfo(CH, PhFileConst.GCI_IMAGECOUNT, pImCount);
    PhGetCineInfo(CH, PhFileConst.GCI_EXPOSURENS, pExpos);
    Exposure=double(pExpos.Value/10^9); % Put Exposure in S
    TimeVector=[0;0;0];
    
    %Subtract the trigger time ( trigger is time zero, not begining of
    %frame one )
    Trigger= libstruct('tagTIME64');
    pTrigger=libpointer('tagTIME64',Trigger);
    pTrigger = libpointer('uint32',0); % initialize a pointer to zero
    PhGetCineInfo(CH, PhFileConst.GCI_TRIGTIMESEC, pTrigger);
    Trigger = pTrigger.Value; % Time in seconds since Epoch
    PhGetCineInfo(CH, PhFileConst.GCI_TRIGTIMEFR, pTrigger);
    Trigger = Trigger + double(int64((pTrigger.Value*10.^8)/2.^32))/10^8;
    % Fractions of a second stored in TRASH units. 
    % The *10^8 is necessary to prevent the /2^32 from zeroing everthing
    % Stores in seconds.
    
    % Loop through the non-negative frames, convert from 16bit to normal
    % Because sometimes the phantomstalker or PCC doesn't actually trim the
    % shot to the save range, we start at firstIm.
    if firstIm<0;firstIm=0;end % This happens sometimes.
    
    %% Loop over frames, save time and pixels
    for frame = firstIm:lastIm
        if mod(frame,100) == 0
            disp(['Frames Index Mod 100: ' num2str(frame) '/' num2str(lastIm)]);
        end
        sFrame.First=frame;
        % Get frame. The 500xxx is the storage buffer size
        [HRES, pixels, IH] = PhGetCineImage(CH,sFrame,(50000000)); 
        % Convert from two eight bit blocks into a single value
        % the second block is the "255'ths" place digit
        pVec = double(pixels(1:2:indNum-1))+225.*double(pixels(2:2:indNum)); 
        CineArray(:,:,frame-firstIm+1) = reshape(pVec,IH.biWidth,IH.biHeight)'+0; %Sometimes there's an offset of 1-3 counts
        %[CineArray(:, :, frame+1), origIm] = PhGetCineImage(CH, frame+uint32(firstIm), false);
        ImageNumber=frame;
        
        % Get time
        [ HRES ] = PhGetCineAuxData( CH, ImageNumber, SelectionCode, DataSize, pEndofFrame);
        
        EndofFrame=pEndofFrame.Value;
        FrameTime = double(int64((EndofFrame.fractions*10.^8)/2.^32))/10^8;
        FrameTime = FrameTime + EndofFrame.seconds; % Seconds just stored as epoch seconds
        TimeVector(frame-firstIm+1)=FrameTime-Trigger; %the 2^32 is because Phantom is stupid (check the documentation for TIME64)
        
        if TimeVector(frame-firstIm+1)/100 >= MovieStop
            break
        end
            
    end
  

    %Subtract half the exposure time ( put stamp in the middle of the
    %exposure, not at end )
    TimeVector=TimeVector-(Exposure/2);
    TimeVector=TimeVector*1e6; % Save in uS
    
    %%  SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %% Save the Cine Matrix and time vector in new folder
    disp('Saving Movie');
    % Trim movie
    while TimeVector ~= sort(TimeVector)
        TimeVector = TimeVector(1:end-1);
    end
    trimInd = find(TimeVector < MovieStop, 1, 'last');
    TimeVector = TimeVector(1:trimInd);
    CineArray = CineArray(:,:,1:trimInd);
     MatFile=sprintf([SavePath 'Shot %d.mat'],shotnum);
%    SaveFile=sprintf('G:\\Cines\\Shot %s.mat',ShotNameVector{j});
     save(MatFile,'CineArray','TimeVector','-v7.3'); %-v7.3 required if file is real big
     
     % Plot the time vector if old data, to verify that its working
     fig=figure; 
     if OldData
         pos=get(fig,'position');
         set(fig,'position',[pos(1:2),pos(3)*2,pos(4)])
         subplot(1,2,2);
         plot(TimeVector);xlabel('Time Point');ylabel('Time [\mu{s}]');
         title('Time Vector');
         subplot(1,2,1);
     end
     surf(sum(CineArray,3)./size(CineArray,3)); shading interp; view([ 0 -90]); colorbar;
     title(['Shot: ' num2str(shotnum)]);
     
     drawnow; % pull buffer immediately.
     %% Save Cine
     if ~OldData
         
         
         %CineFile=sprintf([SavePath '%d.cine'],shotnum);
         CineFile=sprintf([SaveCinePath '%d.cine'],shotnum);
         pFile = libpointer('voidPtr',[int8(CineFile) 0]); % Because C. That's why.
         PhSetCineInfo(CH,PhFileConst.GCI_SAVEFILENAME,pFile);
         sRange=libstruct('tagIMRANGE');
         sRange.Cnt = trimInd;
         sRange.First = 0;
         PhSetCineInfo(CH,PhFileConst.GCI_SAVERANGE,sRange);
         calllib('phfile','PhWriteCineFile',CH,[]); % sychronously save file
         %% Reset Camera
        disp('Save Successful, Resetting Camera');
        PhDestroyCine(CH);
        PhRecordCine(0); % Back to waiting for trigger
        [HRES]=calllib('phcon', 'PhBlackReferenceCI', 0,[]); % Current Session Reference
     end

    %% Save to Tree
    if ~pdc3 && ~OldData
        HitTree = Tree('hitsi3',shotnum,'EDIT');
        HitTree.setDefault(HitTree.getNode('\HITSI3::TOP.IMAGING3'))
        try
            HitTree.deleteNode('FASTCAM');
        end
        HitTree.addNode('FASTCAM','NUMERIC');
        curr=HitTree.getNode('FASTCAM');
        curr.addTag('\FASTCAM');
        for i = 1:length(TimeVector)
            currTime = TimeVector(i);
            FcurrTime = Float32(currTime);
            seg = Float32Array(CineArray(:,:,i));
            curr.makeSegment(FcurrTime,FcurrTime,Float32Array(currTime),seg);
        end
        HitTree.write;
    end
    
%{
    %% Saving all the stuff we want into the tree
    MovieStopX = find(TimeVect or>0 & TimeVector<MovieStop); %this trims the cine array
    XArray = CineArray(:,:,MovieStopX);
    nFrames = numel(MovieStopX);
    nTime = numel(TimeVector);
    interval = max(TimeVector)/nTime;%this is dt b/w frames

    mdsconnect('landau.hit')
    mdsopen('imaging3', 0);
    mdsvalue('BeginSegment(\phantom,$1,$2,makedim(*:$1:$2:$3))',0,MovieStop,interval,uint16(XArray));

     for n = 1:nFrames;
         mdsvalue('PutSegment(\phantom, -1,$)',uint16(XArray(:,:,n)));
     end

    mdsput('\imaging3_Image_Height', Height);
    mdsput('\imaging3_Image_Width', Width);
    mdsput('\imaging3_Image_Exposure', ExposureTime);
    mdsclose;
    mdsdisconnect;
%}

end







% NOTES:
% check Phantom documentation for CINE specific functions ( its included in 
% the file path for this program on the Cine computer
% Get the Exposure time from the Phantom 692 program on the cine computer
% To test the cine array, on matlab 2013, implay(~CineArray);
% I dont know why we cant get the exposure normally
% Code works on R11 and R13, maynot work below R7 ( phantom compatability issues
% may arrise at that point.

  