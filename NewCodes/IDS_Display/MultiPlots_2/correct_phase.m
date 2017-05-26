%% Correct phase from sine_fit and make sanity plots, for Multiplots_2
function [dataPhase] = correct_phase(dataPhase,plotSanityPhase,n,i,ax,in,param,data,impacts)
    % mod the phase
    % Sanity check
    if plotSanityPhase==1;eval(['plot(ax(' num2str(9+n) '),dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth)']);end
    dataPhase(:,i) = mod(dataPhase(:,i),2*pi);
    if plotSanityPhase==1;eval(['plot(ax(' num2str(9+n) '),dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'','':'')']);end
    % test automatic phase finding
    for j = 1: size(param,1) % for all impacts
        dataPhase(j,i) = find_Phase(param(j,1+(5*(i-1)): 5 + 5*(i-1),n));
    end
    if plotSanityPhase==1;eval(['plot(ax(' num2str(9+n) '),dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''--'')']);end
   % shift lines if a periodicity jump occurs
    for j = 2:size(data, 2)
        if (dataPhase(j,i) - dataPhase(j-1,i))>(pi)
            dataPhase(j,i) = dataPhase(j,i)-(2*pi);
            disp(['JUMP: -2Pi, Line: ' num2str(i) ' Impact: ' num2str(impacts(j))]);
        elseif (dataPhase(j,i) - dataPhase(j-1,i))<(-pi)
           disp(['JUMP: +2Pi, Line: ' num2str(i) ' Impact: ' num2str(impacts(j))]);
            dataPhase(j,i) = dataPhase(j,i)+(2*pi);
        end
    end
    if plotSanityPhase==1;eval(['plot(ax(' num2str(9+n) '),dat(1).impacts(1:size(data,2)),dataPhase(:,i).*180./pi,''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''-.'')']);end
    if plotSanityPhase==1;
        legend(ax(10),{'Raw Phase';'Mod 2\pi Phase';'Find\_Phase';'Jump Corrected';},'location','SouthWest');
        set(ax(11),'xticklabel',[]);
        set(ax(12),'xticklabel',[]);
        eval(['ylabel(ax' num2str(9+n) ','' Phase [deg]'')']);
        xlabel(ax(10),'Impacts [cm]');
        mBox1=uicontrol(h(10),'Style','Text');% Make the title
        set(mBox1,'String',['Phase Calculation Methods']); 
        set(mBox1,'fontweight','bold'); set(mBox1,'fontsize',13);
        set(mBox1,'Position', [ 80,960,450,20]);set(mBox1,'BackgroundColor',[1 1 1]);
        eval(['title(ax(' num2str(9+n) '), ''Shot: ' num2str(in(n).shot) ''')']);
    end
    dataPhase(:,i) = mod(param(:,4+5*(i-1),n),2*pi); % 2Pi Periodicity
    
end