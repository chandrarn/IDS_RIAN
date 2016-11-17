%%%CINE-TO-MAT CONVERTER PROGRAM
%%RIAN CHANDRA 2014 & KYLE ROBERTS 2016 & RIAN CHANDRA 2016
% open cine file from camera
%         read in cine
%         make time vector
%         Save .mat file
%         Throw to tree
%         Return Camera to Initial State
 
%% The idea is to have some kind of while loop that allows the program to run in the background and wait for a shot.
clear all;
%close all;
%clc;
import MDSplus.*
addpath(genpath('T:\RChandra\Phantom\PhMatlabSDK'));
pdc3 =0;
OldData = 1;
shots =  [160412018:160412022];[160518015,160518017,160518027,160518029];% must be in ascending order
SavePath = 'T:\\PhantomMovies\\';
LoadPhantomLibraries();
RegisterPhantom(true);
%% Getting time signatures associated with begining and ending of movie
if ~OldData
    Conn = Connection('landau.hit');
    Conn.openTree('hitsi3',0);
    MovieStop = double(Conn.get('\Shot_Length')) + 0.0005;
    MovieStop = MovieStop*10^(6);
    HRES=calllib('phcon', 'PhBlackReferenceCI', 0,[]); % CSR
    if HRES == 132
        disp('Successfully refreshed CCD');
    elseif HRES == -259
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
    MovieStop = 10*10^6;
end
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
cont = 1;shotnum=0;
emergencycounter = 0;

while cont==1;
    display('Waiting for store to conclude');
    % Wait for the shot to happen
    if ~OldData
        Conn.openTree('hitsi3',0); % will triger for pdc3 and hitsi3
        Conn.get('wfevent("END_OF_STORE")');
    end
    if OldData
            shotnum=shots(find(shots>=(shotnum+1),1,'first'));
            if shotnum==shots(end) % if we've reached the end, stop
                cont=0;
            end
    else
        try
            Conn.openTree('pdc3',0);
            pdcNew = double(Conn.get('$Shot'));
        catch 
            pdcNew = NaN;
        end
        Conn.openTree('hitsi3',0);
        hitsiNew = double(Conn.get('$Shot'));
        if hitsiNew ~= hitsiShot
            shotnum = hitsiNew;
            pdc3=0; % set flag for storing to tree;
        elseif ~(isnan(pdcNew)) && (pdcNew ~= pdc3Shot) 
            shotnum = pdcNew;
            pdc3=1;
        else
            disp('Failed to get new shot number, trying again');
            shotnum = hitsiShot + emergencycounter;
            emergencycounter=emergencycounter+1;
            %continue; % Skip rest of program, wait for next shot
        end
        hitsiShot=hitsiNew; % update the shot number
        pdc3Shot = pdcNew;
    end
    
    %clearvars -except 'Exposure' 'Shots' 'j' 'ShotNameVector' 'date' ;
    % Figure out how to get the Shot number
    fprintf('Working on Cine #%d...\n', shotnum)
    
    %This creates a new handle for the cine 1, in camera 0. The handle is CH
    if ~OldData
        [HRES, CH] = PhNewCineFromCamera(0,1);
    else
        [HRES, CH] = PhNewCineFromFile(['T:\PhantomMovies\' num2str(shotnum) '.cine']);
    end
    
    % Reset Image Processing Parameters
    PhSetCineInfo(CH,PhFileConst.GCI_BRIGHT,libpointer('float',0.0));
    PhSetCineInfo(CH,PhFileConst.GCI_CONTRAST,libpointer('float',1.0));
    PhSetCineInfo(CH,PhFileConst.GCI_GAMMA,libpointer('float',1.0));
    
%%  %% Put Cine into matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    indNum = IH.biSizeImage;
    
    
%     %Make time vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    TimeVector=[0;0;0];
    
    %Subtract the trigger time ( trigger is time zero, not begining of
    %frame one )
    Trigger= libstruct('tagTIME64');
    pTrigger=libpointer('tagTIME64',Trigger);
    PhGetCineInfo(CH, PhFileConst.GCI_TRIGTIMEFR, pTrigger);
    Trigger=pTrigger.Value;   
    Trigger = double(int64((Trigger.fractions*10.^8)/2.^32)); % preform subtraction
    % Loop through the non-negative frames, convert from 16bit to normal
    for frame = 0:lastIm
        if mod(frame,100) == 0
            disp(['Frames Index Mod 100: ' num2str(frame) '/' num2str(lastIm)]);
        end
        sFrame.First=frame;
        [HRES, pixels, IH] = PhGetCineImage(CH,sFrame,(50000000));
        % Convert from two eight bit blocks into a single value
        pVec = double(pixels(1:2:indNum-1))+225.*double(pixels(2:2:indNum)); 
        CineArray(:,:,frame+1) = reshape(pVec,IH.biWidth,IH.biHeight)'+0; %Sometimes there's an offset of 1-3 counts
        %[CineArray(:, :, frame+1), origIm] = PhGetCineImage(CH, frame+uint32(firstIm), false);
        ImageNumber=frame;
        [ HRES ] = PhGetCineAuxData( CH, ImageNumber, SelectionCode, DataSize, pEndofFrame);
        
        EndofFrame=pEndofFrame.Value;
        TimeVector(frame+1)=double(int64((EndofFrame.fractions*10.^8)/2.^32))-Trigger; %the 2^32 is because Phantom is stupid (check the documentation for TIME64)
        
        if TimeVector(frame+1)/100 >= MovieStop
            break
        end
            
    end
    
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Make time vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp('Creating Timebase');
%     %Create variables to select and store the end-of-frame times
%     EndofFrame = libstruct('tagTIME64');
%     pEndofFrame = libpointer('tagTIME64', EndofFrame);
%     SelectionCode = PhConConst.GAD_TIME;
% %     ExposureTime = PhConConst.GAD_EXPOSURE;
%     DataSize = uint32(32);
%     exposure = uint16(16);
%     %Get total number of frames
%     pImCount = libpointer('uint32Ptr',0);
%     pExpos = libpointer('uint32Ptr',0);
%     PhGetCineInfo(CH, PhFileConst.GCI_IMAGECOUNT, pImCount);
%     PhGetCineInfo(CH, PhFileConst.GCI_EXPOSURENS, pExpos);
%     TimeVector=[0;0;0];
%     for index = 0:lastIm-1 % pull the end of frame times for each frame
%         if mod(frame,100) == 0
%         disp(['TimeBase Index Mod 100: ' num2str(index)]);
%         end
%         ImageNumber=index;
%         [ HRES ] = PhGetCineAuxData( CH, ImageNumber, SelectionCode, DataSize, pEndofFrame);
%         
%         EndofFrame=pEndofFrame.Value;
%         TimeVector(index+1)=double(int64((EndofFrame.fractions*10.^8)/2.^32)); %the 2^32 is because Phantom is stupid (check the documentation for TIME64)
%    
%     end
%    
%         %Subtract the trigger time ( trigger is time zero, not begining of
%     %frame one )
%     Trigger= libstruct('tagTIME64');
%     pTrigger=libpointer('tagTIME64',Trigger);
%     PhGetCineInfo(CH, PhFileConst.GCI_TRIGTIMEFR, pTrigger);
%     Trigger=pTrigger.Value;   
%     TimeVector=TimeVector-double(int64((Trigger.fractions*10.^8)/2.^32)); % preform subtraction
    %Subtract half the exposure time ( put stamp in the middle of the
    %exposure, not at end )
    Exposure=double(pExpos.Value/2/1000); % Put Exposure in uS
    TimeVector=TimeVector-Exposure;
    TimeVector=TimeVector/100;
    Exposure = Exposure/50;
    
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
     figure; surf(sum(CineArray,3)./size(CineArray,3)); shading interp; view([ 0 90]); colorbar;
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

  