function out = calSecLine % find the offset of the second fiber
    out(:,2) = [53.5,-1.4/100,.001,-.1];
    out(:,3) = [52.5,-2/100,-.1,-.5]'
    out(:,4) = [54,-.5/100,.1,.5];
    
    out(:,1) = lsqcurvefit(@SecDegreeFit,out(:,2),dat(1).param.peaks(:,2),dat(1).param.peaks(:,3),out(:,3),out(:,4));
    %plot data
    figure; plot(dat(1).param.peaks(:,2),dat(1).param.peaks(:,3),'*',dat(1).param.peaks(:,2),SecDegreeFit(out(:,1),dat(1).param.peaks(:,2)))
    title(['2^{nd} Line Scalar Offset: ' num2str(out(4,1))]);
    legend({'Raw Data';'Linear Fit'});
end
function y = SecDegreeFit(par,x)
% par(1) x^0 coeffecient
% par(2) x^1 coeffecient
% par(3) x^2 coeffecient
% par(4) 2nd bundle offset (29:end)

temp = x;
%temp(29:end) = temp(29:end)+par(4); % apply offset to second bundle 

y = par(1) + par(2).*(temp) + par(3).*(temp.^2); % value of the line at every point.

y(30:end) = y(30:end) +par(4); % apply offset to second set of fibers.

end
