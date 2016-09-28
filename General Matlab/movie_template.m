%% call the movie maker
clear all
close all
%Set movie properties
vidName = 'movie.avi';
framerate = 5;
%imgname_base is the root of your image files for the frames
imgname_base = 'S:\Nimrod\nimfl\Movie\hi-res\nimfl0';
%create video object
video = VideoWriter(vidName,'Motion JPEG AVI');
video.FrameRate = framerate;
video.Quality = 75;
open(video);

%frame increment
dumpstart = 100;
dumstep = 50;
dumpend = 1000;

%put preparations for making the plots here


for it = 1:length(dumps)
    %create image file name for the frame
    %imgname_base with the frame number appended to the end
    fpng = [imgname_base,int2str(dumps(it)),'.png'];
    try %to read the image of the frame if its already been made
        g = imread(fpng,'png');
    catch %try writing the image file
        try %to create the plot
            f1 = figure('visible','off');
            %Put plotting routine here
            
            %Save file as a png
            print(f1,'-dpng',fpng);
            %delete figure object so matlab doesn't get confused
            delete(f1);
            %read image file
            g = imread(fpng,'png');
        catch %creation of plot failed
            fprintf('Unable to create frame for image %i\n',dumps(it));
            %commenting out continue just makes a new frame with the
            %previous frame, otherwise it skips to try the next frame
            continue
        end
    end
    %write frame to image file
    writeVideo(video,g);
    fprintf('Frame for image %i added.\n',dumps(it));
end
%End by closing the video file
close(video);