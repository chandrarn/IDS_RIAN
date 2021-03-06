function out = calSecLine1(data,index) % find the offset of the second fiber
    % Allows second fiber to have different slope. Which happens,
    % sometimes?
    
    %PEAKS = dat(1).param.peaks; %for old data
    PEAKS = data; % for calMain
    
    f=@(par,x)SecDegreeFit(par,x,index);
    
    % automatically find intial guess and bounds
    offset = PEAKS(1,3)+.5;
    slope1 = (PEAKS(25,3)-PEAKS(5,3))/(PEAKS(25,2)-PEAKS(5,2));
    slope2 = (PEAKS(65,3)-PEAKS(45,3))/(PEAKS(65,2)-PEAKS(45,2));
    shift = PEAKS(45,3)-(offset+slope2*PEAKS(45,2));
    out(:,2) = [offset,slope1,slope2,.001,shift];
    out(:,3) = [offset-1.5,[slope1,slope2]*10,-.1,shift-.5];
    out(:,4) = [offset+1.5,[slope1,slope2]*.1,.1,shift+.5];
    if slope1>0
        out(2:3,3)=out(2:3,3)./100;
        out(2:3,4)=out(2:3,4).*100;
    end
    
    %out(5,:) = index; % where to apply the shift

    % valid for 160518,19
%     out(:,2) = [82,-2/100,.001,-.1];
%     out(:,3) = [81.5,-3/100,-.1,-.5]; %lb
%     out(:,4) = [82.5,-.5/100,.1,.5]; % ub
    out(:,1) = lsqcurvefit(f,out(:,2),PEAKS(:,2),PEAKS(:,3),out(:,3),out(:,4));

    %plot data
    figure; plot(PEAKS(:,2),PEAKS(:,3),'*',PEAKS(:,2),SecDegreeFit(out(:,1),PEAKS(:,2),index))
    title(['2^{nd} Line Scalar Offset: ' num2str(out(5,1))]);
    legend({'Raw Data';'Linear Fit'});
    if max(ismember(out(:,1),out(:,3:4)))
        assignin('base','out',out);
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

y(1:index-1) = par(1) + par(2).*(x(1:index-1)) + par(4).*(x(1:index-1).^2); % value of the line at every point.
y(index:length(x)) = par(1) + par(3).*(x(index:end)) + par(4).*(x(index:end).^2);
try
    y(index:end,1) = y(index:end,1) +par(5); % apply offset to second set of fibers.
catch % sometimes this fails for really unclear reasons, 30 is usually where the break is: 151217026
    y(43:end,1) = y(43:end,1) +par(5);
end
y = y';
end
