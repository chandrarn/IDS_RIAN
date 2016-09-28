%
%===================================================================
%
%Adam Madson
%1/05/05
%
%ABEL INVERSION CALL FUNCTION: abel1.m
%
%NOTE: Requires centroid.m and centroiderror.m to run.  Also uses
%   magneticmode.m data for some plots.
%
%Description: The Abel routine determines the radial plasma profile
%   from the chord integrated data using the Abel inversion method.
%   The abel1.m and abel2.m routines are identical except that they
%   use PDA1 or PDA2 data.  Movie and plotting routines are listed
%   at the end of this code.
%
%===================================================================
%
% begin function
function [Lshell1,Rshell1] = abel1(X1,nonzero1,T1t,Centroid1t,...
       PDA1t,...                      %Chord data used
       magtime,mode1,Cx,Cy,...        %Mode data for subplot
       PDA1,shotnum,Tms_1,...         %Contour data for subplot
       PDA1use,m1,moviemake);                   %PDA in use

%function without magnetic data
% function [Lshell1,Rshell1] = abel1(X1,nonzero1,T1t,Centroid1t,...
%        PDA1t,...                      %Chord data used
%        PDA1,shotnum,Tms_1,...         %Contour data for subplot
%        PDA1use);                      %PDA in use
       
format short

%BREAK IF PDA NOT IN USE
for alpha = 1
    if PDA1use  == 0
        Lshell1=0; Rshell1=0;
        break
    else    

%vector and array creation
n = length(T1t);
leftindex = [];
rightindex = [];
left = [];
right = [];
Lvect = [];
Rvect = [];
dleft = [];
dright = [];
Lshell1 = zeros(length(nonzero1),n);
Rshell1 = zeros(length(nonzero1),n);

%Ascending order check
if X1(1) < X1(2)
    X1 = X1;
    PDA1t = PDA1t;
elseif X1(1) > X1(2)
    X1 = fliplr(X1);
    PDA1t = flipud(PDA1t);
end
    
%create length vector to define shells by distance from centroid
for i = 1:n
    signal = PDA1t(:,i);
    
    %left side
    leftindex = find(X1 <= Centroid1t(i));
    left = X1(leftindex);
    Lvect = abs(left-Centroid1t(i));%chord distances from centroid
    for p = 1:length(Lvect)-1
        dleft(p) = abs(Lvect(p+1)-Lvect(p));
    end
    dleftavg = mean(dleft);
    %shell radii from centroid and 1st chord distance from cent.
    leftvect = [Lvect(1)+dleftavg, Lvect];  
                                            
    leftsignal = signal(leftindex);         %signal at each chord
    L1 = zeros(length(leftvect)-1);
    L2 = zeros(length(leftvect)-1);
    
    for j = 1:length(leftvect)-1
        for k = 1:length(leftvect)-1
            if leftvect(k)^2-leftvect(j+1)^2 > 0
                L1(j,k) = sqrt(leftvect(k)^2-leftvect(j+1)^2);
            else
                L1(j,k) = 0;
            end
        end
        for k = 1:length(leftvect)-1
            if leftvect(k)^2-leftvect(j)^2 > 0
                L2(j,k) = sqrt(leftvect(k+1)^2-leftvect(j+1)^2);
            else
                L2(j,k) = 0;
            end
        end
        L = 2*(L1-L2);                  %left side length matrix
    end
    LS = L\leftsignal;
    
    %Left side shell emissivities at given time steps.  Row order 
    %   from outside to nearest to centroid
    Lshell1(1:length(LS),i) = LS;          
                                        
    %right side
    rightindex = find(X1 >= Centroid1t(i));  
    right = X1(rightindex);
    Rvect = abs(right-Centroid1t(i));%chord distances from centroid
    for p = 1:length(Rvect)-1
        dright(p) = abs(Rvect(p+1)-Rvect(p));
    end
    drightavg = mean(dright);
    %shell radii from centroid and 1st chord distance from cent.
    rightvect = [Rvect, Rvect(end)+drightavg]; 
 
    rightsignal = signal(rightindex);       %signal at each chord
    R1 = zeros(length(rightvect)-1);
    R2 = zeros(length(rightvect)-1);
    
    for j = 1:length(rightvect)-1
        for k = 1:length(rightvect)-1
            if rightvect(k+1)^2-rightvect(j)^2 > 0
                R1(j,k) = sqrt(rightvect(k+1)^2-rightvect(j)^2);
            else
                R1(j,k) = 0;
            end
        end
        for k = 1:length(rightvect)-1
            if rightvect(k)^2-rightvect(j)^2 > 0
                R2(j,k) = sqrt(rightvect(k)^2-rightvect(j)^2);
            else
                R2(j,k) = 0;
            end
          end
        R = 2*(R1-R2);          %right side length matrix
    end
    RS = R\rightsignal;
    
    %Right side shell emissivities at given time steps.  Row order 
    %   from outside to nearest to centroid
    Rshell1(length(nonzero1)-length(RS)+1:end,i) = RS;

end

%===================================================================
%%% MOVIE %%%
%===================================================================
% 
%Contour Plot Truncation
% if moviemake
%     first = find(Tms_1>19 & Tms_1<20);
%     last = find(Tms_1>100 & Tms_1<102);
%     [A1,B1] = meshgrid( Tms_1(first(1):last(1)), X1 );
% 
%     %Truncation of Lshell and Rshell to sample every 5th location
%     Lshellmovie = zeros(m1,floor(length(Lshell1)/5));
%     Rshellmovie = zeros(m1,floor(length(Rshell1)/5));
%     Tmovie = zeros(1,5/8*floor(length(T1t)/5));
%     PDA_outmovie = zeros(m1,floor(length(PDA1t)/5));
%     for i = 1:1:length(Tmovie)
%         Lshellmovie(:,i) = Lshell1(:,i*5);
%         Rshellmovie(:,i) = Rshell1(:,i*5);
%         Tmovie(i) = T1t(:,i*5);
%         PDA_outmovie(:,i) = PDA1t(:,i*5);
%     end
% 
%     %mov = avifile('PDA_Abel1.avi','compression','Cinepak')
%     % mov = avifile('test.avi','compression','Cinepak')
%     mov = avifile(['PDA1_Abel',num2str(shotnum),'.avi'],'compression','Cinepak')
%     for i = 1:length(Tmovie)
%         XLindex = find(Lshellmovie(:,i));
%         XRindex = find(Rshellmovie(:,i));
%         XL = X1(XLindex);
%         XR = X1(XRindex);
%         Lshellvalue = Lshellmovie(XLindex,i);
%         Rshellvalue = Rshellmovie(XRindex,i);
% 
%         figure(6)
%         subplot(2,2,1)
%         pcolor(A1, B1, PDA1(:, first(1):last(1)) );shading interp;
%         axis([20 100 -2.5 2.5]);
%         xlabel('Time (\mus)')
%         ylabel('X-Axis Location (cm)')
%         title(['Emissivity Contour for Shot #',int2str(shotnum) ])
%         hold on
%         plot(T1t,Centroid1t,'.','MarkerEdgeColor','k','MarkerSize',4.5)
%         plot([Tmovie(i),Tmovie(i)],[-2.5,2.5],'w');
%         hold off
% 
%         subplot(2,2,2);
%         plot(XL,Lshellvalue,'-.r*',XR,Rshellvalue,'-.b*',...
%             X1,PDA_outmovie(:,i),'k');
%         text(1,.45,['T = ' num2str(Tmovie(i))]);
%         axis([-2.5 2.5 0 1]);
%         xlabel('X-Axis Location (cm)')
%         ylabel('Attenuated Signal')
%         title(['Abel Inversion for Shot #',int2str(shotnum) ])
% 
%         subplot(2,2,3:4);
%         plot(T1t,Centroid1t,'.','MarkerEdgeColor','k','MarkerSize',4.6);
%         hold on
%         plot(magtime,Cx,'r','LineWidth',2);
%         plot([Tmovie(i),Tmovie(i)],[-2.5,2.5]);
%         hold off
%         axis([20 70 -2.5 2.5]);
%         legend('PDA Centroid','Magnetic Mode Centroid','Location','NorthEast')
%         xlabel('Time (\mus)')
%         ylabel('Array Location (cm)')
% 
%         F = getframe(gcf);
%         mov = addframe(mov,F);
%     end
%     mov = close(mov);
% end
% %===================================================================
% %%% DATA PLOTS %%%
% %===================================================================
% %THESIS DATA PLOTS 
% figure(17)
%     
% subplot(2,2,1);           %data at 25 microseconds
% index = find(T1t>25);
% XLindex = find(Lshell1(:,index(1)));
% XRindex = find(Rshell1(:,index(1)));
% XL = X1(XLindex);
% XR = X1(XRindex);
% Lshellvalue = Lshell1(XLindex,index(1));
% Rshellvalue = Rshell1(XRindex,index(1));
% plot(XL,Lshellvalue,'-.r*',XR,Rshellvalue,'-.b*',X1,...
%       PDA1t(:,index(1)),'-k*');
% title(['Time = ' num2str(T1t(index(1))) '\mus'])
% xlabel('X-Axis Location (cm)')
% ylabel('Attenuated Signal')
% axis([-2.5 2.5 0 1]);
% 
% subplot(2,2,2);           %data at 45 microseconds
% index = find(T1t>40);
% XLindex = find(Lshell1(:,index(1)));
% XRindex = find(Rshell1(:,index(1)));
% XL = X1(XLindex);
% XR = X1(XRindex);
% Lshellvalue = Lshell1(XLindex,index(1));
% Rshellvalue = Rshell1(XRindex,index(1));
% plot(XL,Lshellvalue,'-.r*',XR,Rshellvalue,'-.b*',X1,...
%       PDA1t(:,index(1)),'-k*');
% title(['Time = ' num2str(T1t(index(1))) '\mus'])
% xlabel('X-Axis Location (cm)')
% ylabel('Attenuated Signal')
% axis([-2.5 2.5 0 1]);
% 
% subplot(2,2,3);           %data at 65 microseconds
% index = find(T1t>60);
% XLindex = find(Lshell1(:,index(1)));
% XRindex = find(Rshell1(:,index(1)));
% XL = X1(XLindex);
% XR = X1(XRindex);
% Lshellvalue = Lshell1(XLindex,index(1));
% Rshellvalue = Rshell1(XRindex,index(1));
% plot(XL,Lshellvalue,'-.r*',XR,Rshellvalue,'-.b*',X1,...
%       PDA1t(:,index(1)),'-k*');
% title(['Time = ' num2str(T1t(index(1))) '\mus'])
% xlabel('X-Axis Location (cm)')
% ylabel('Attenuated Signal')
% axis([-2.5 2.5 0 1]);
% 
% subplot(2,2,4);           %data at 85 microseconds
% index = find(T1t>60);
% XLindex = find(Lshell1(:,index(end)));
% XRindex = find(Rshell1(:,index(end)));
% XL = X1(XLindex);
% XR = X1(XRindex);
% Lshellvalue = Lshell1(XLindex,index(end));
% Rshellvalue = Rshell1(XRindex,index(end));
% plot(XL,Lshellvalue,'-.r*',XR,Rshellvalue,'-.b*',X1,...
%       PDA1t(:,index(end)),'-k*');
% title(['Time = ' num2str(T1t(index(end))) '\mus'])
% xlabel('X-Axis Location (cm)')
% ylabel('Attenuated Signal')
% axis([-2.5 2.5 0 1]);
%
% %SAVE TO .jpeg FILE
% filename = (['L:\Users\Adam\PDA\Collected Matlab Data\',num2str(shotnum)...
%     ' Abel1.jpeg']);
% print('-f7', '-djpeg', filename);
% 
% SAVE TO .eps FILE FOR HARD COPY PLOTS 
% filename = (['L:\Users\Adam\PDA\Collected Matlab Data\',num2str(shotnum)...
%    ' Abel1.eps']);
% print('-f7', '-depsc', filename);

%PDA USE LOOP
    end
end

