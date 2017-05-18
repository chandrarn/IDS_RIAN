%% Plot Data Reconstruction with Fits for Multiplots_2
function [] = plotReconst(dat,in,param,n,i,ax,h,plotChan)

if in(1).shot >129500
    %if i==2 % only plot the second time around
    %plotChan=9;
    for j = 1:2 % only plot the special explainatory reconstruction once
        reconTime = linspace(dat(1).time(1).*in(n).timeScale,dat(1).time(end).*in(n).timeScale,length(dat(1).time)*10);
        eval(['plot(ax(' num2str(12+n+(j==2)*1) '),dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan),''color'',partialColor{n,1},''linewidth'',lnwdth)']);
        eval(['plot(ax(' num2str(12+n+(j==2)*1) '),dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan+(length(dat(1).impacts))/2),''color'',partialColor{n,2},''linewidth'',lnwdth,''linestyle'',''-'')']);
        eval(['plot(ax(' num2str(12+n+(j==2)*1) '),reconTime,param(plotChan,2+5*(i-1),n)+param(plotChan,3+5*(i-1),n)*sin(2*pi*reconTime.*1e-3.*14500 + dataPhase(plotChan,i)),''color'',partialColor{n,i},''linewidth'',lnwdth,''linestyle'',''--'')']);
        %eval(['plot(ax' num2str(12+n) ',dat(1).time.*in(n).timeScale,param(plotChan,7,n)+param(plotChan,8,n)*sin(2*pi*dat(1).time.*1e-3.*in(n).timeScale.*14500 + dataPhase(plotChan,2)),''color'',partialColor{n,2},''linewidth'',lnwdth,''linestyle'',''--'')']);
        eval(['title(ax9' num2str(12+n) '), ''Shot: ' num2str(in(n).shot) ''')']);
        set(ax(15),'xticklabel',[]);
        set(ax(14),'xticklabel',[]);
        eval(['ylabel(ax(' num2str(12+n+(j==2)*1) '),'' Velocity [km/s]'')']);
        xlabel(ax(13),'Time [ms]');
        mBox2=uicontrol(h(11),'Style','Text');% Make the title
        set(mBox2,'String',['Raw Data and Reconstruction, Impact: ' num2str(dat(1).impacts(plotChan)) '[cm]']); 
        set(mBox2,'fontweight','bold'); set(mBox2,'fontsize',13);
        set(mBox2,'Position', [ 80,960,450,20]);set(mBox2,'BackgroundColor',[1 1 1]);
        % Explainatory plot stuff
        if j==2
            plot(ax(16),[reconTime(1)-.025,reconTime(end)+.025],ones(1,2).*param(plotChan,2+5*(i-1),n),'k--')
            plot(ax(16),ones(1,2)*(28/14500 +.25/14500 -dataPhase(plotChan,i)/(2*pi*14500))*1e3,[1,13],'k--')
            title(ax(16),['Raw Data and Fit, Shot: ' num2str(in(n).shot) ', Impact: ' num2str(dat(1).impacts(plotChan)) ' [cm]']);
            %legend(ax16,{['Upper Fiber, RMS Error: ' num2str(RMS(1,plotChan,3,10)) ' [km/s]'];['Lower Fiber, RMS Error: ' num2str(RMS(2,plotChan,3)) ' [km/s]']},'Location','SouthWest','fontsize',10)
            xlabel(ax(16),'Time [ms]');
        end


    end

    % Explainatory reconsturction
else
    plotChan=14;
    eval(['plot(ax(' num2str(12+n) '),dat(1).time.*in(n).timeScale,dat(in(n).line).vel(:,plotChan),''color'',partialColor{n,1},''linewidth'',lnwdth)']);
    eval(['plot(ax(' num2str(12+n) '),dat(1).time.*in(n).timeScale,param(plotChan,2,n)+param(plotChan,3,n)*sin(2*pi*dat(1).time.*1e-3.*in(n).timeScale.*14500 + dataPhase(plotChan)),''color'',partialColor{n,1},''linewidth'',lnwdth,''linestyle'',''--'')']);
    set(ax(15),'xticklabel',[]);
    set(ax(14),'xticklabel',[]);
    eval(['title(ax(' num2str(12+n) '), ''Shot: ' num2str(in(n).shot) ''')']);
    eval(['ylabel(ax(' num2str(12+n) '),'' Velocity [km/s]'')']);
    xlabel(ax(13),'Time [ms]');
    mBox2=uicontrol(h(11),'Style','Text');% Make the title
    set(mBox2,'String',['Raw Data and Reconstruction, Impact: ' num2str(dat(1).impacts(plotChan)) '[cm]']); 
    set(mBox2,'fontweight','bold'); set(mBox2,'fontsize',13);
    set(mBox2,'Position', [ 80,960,450,20]);set(mBox2,'BackgroundColor',[1 1 1]);
end

linkaxes([ax(10),ax(11),ax(12)],'x');
linkaxes([ax(13),ax(14),ax(15)],'x');

end