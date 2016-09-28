function out = calSecLine(data,index) % find the offset of the second fiber
    
    %PEAKS = dat(1).param.peaks; %for old data
    PEAKS = data; % for calMain
    
    f=@(par,x)SecDegreeFit(par,x,index);
    
    % automatically find intial guess and bounds
    offset = PEAKS(1,3)+.5;
    slope = (PEAKS(25,3)-PEAKS(5,3))/(PEAKS(25,2)-PEAKS(5,2));
    shift = PEAKS(33,3)-(offset+slope*PEAKS(33,2));
    out(:,2) = [offset,slope,.001,shift];
    out(:,3) = [offset-1.5,slope*10,-.1,shift-.5];
    out(:,4) = [offset+1.5,slope*.1,.1,shift+.5];
    
    %out(5,:) = index; % where to apply the shift

    % valid for 160518,19
%     out(:,2) = [82,-2/100,.001,-.1];
%     out(:,3) = [81.5,-3/100,-.1,-.5]; %lb
%     out(:,4) = [82.5,-.5/100,.1,.5]; % ub
    out(:,1) = lsqcurvefit(f,out(:,2),PEAKS(:,2),PEAKS(:,3),out(:,3),out(:,4));

    %plot data
    figure; plot(PEAKS(:,2),PEAKS(:,3),'*',PEAKS(:,2),SecDegreeFit(out(:,1),PEAKS(:,2)))
    title(['2^{nd} Line Scalar Offset: ' num2str(out(4,1))]);
    legend({'Raw Data';'Linear Fit'});
    if max(ismember(out(:,1),out(:,2:4)))
        error('ERROR: calSecLine hit a bound');
    end
end
function y = SecDegreeFit(par,x,index)
% par(1) x^0 coeffecient
% par(2) x^1 coeffecient
% par(3) x^2 coeffecient
% par(4) 2nd bundle offset (30:end)
% index fiber break

%temp = x;
%temp(29:end) = temp(29:end)+par(4); % apply offset to second bundle 

y = par(1) + par(2).*(x) + par(3).*(x.^2); % value of the line at every point.

try
    y(index:end,1) = y(index:end,1) +par(4); % apply offset to second set of fibers.
catch % sometimes this fails for really unclear reasons, 30 is usually where the break is: 151217026
    y(30:end,1) = y(30:end,1) +par(4);
end

end
