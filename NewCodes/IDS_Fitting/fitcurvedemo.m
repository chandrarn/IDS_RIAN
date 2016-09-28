% function [estimates, model] = fitcurvedemo(xdata, ydata)
% % Call fminsearch with a random starting point.
% start_point = [10.0,60.0,3.14];
% model = @expfun;
% estimates = fminsearch(model, start_point,optimset('maxFunEvals',1e+10))
% % expfun accepts curve parameters as inputs, and outputs sse,
% % the sum of squares error for A*exp(-lambda*xdata)-ydata,
% % and the FittedCurve. FMINSEARCH only needs sse, but we want
% % to plot the FittedCurve at the end.
%     function [sse, FittedCurve] = expfun(params)
%         A = params(1);
%         f = params(2);
%         phi=params(3);
%         FittedCurve = A .* sin(2*pi*f*xdata+phi);
%         ErrorVector = FittedCurve - ydata;
%         sse = sum(ErrorVector .^ 2);
%     end
% end
% 



        
function [ bestcoeffs ]= fitcurvedemo(dat)

 close all;
        addpath ('T:\IDS\Data Repository');
        shots = [12858710];
        addpath('T:\RChandra\A-A-Ron Code\General Matlab\extrema');
        TimeBounds= [1,2];
       
        %how to guess frequency: yPos=y.*(y>0); Eliminates times when the
        %graph spkies randomly in the negative region
        %sort imax, disregard first entry and last, take two consecutive
        %values, get x(highervalue)-x(lowervalue) and that should be the
        %period in miliseconds.
        
        
        %Guessing phase: if the value is closer to zero than to A and goin
        %g down, its sifted by pi, if goin gup, by 3/2pi, if closer to A
        %and going down, by 2/3pi, going up, by 1/3pi, if negative A: going
        %down, 4/3 pi, going up, 5/3. I think. Math is hard
        %take final result and mod (2pi) to normalize
        
        
        
        % eval(sprintf('load(''dat%i'');', 12858710)); % get the shot
         TimeWing = 1/20;
         minBound=find(dat.iinjxTime>(1.75-.000001 - TimeWing));
        maxBound=find(dat.iinjxTime>(1.75-.000001 + TimeWing));
        minBound=minBound(1); maxBound=maxBound(1);
        x=linspace(dat.iinjxTime(minBound),dat.iinjxTime(maxBound),501)';
        y=dat.iinjx(minBound:maxBound);
        
        %guess amplitude
        A=max(dat.iinjx(minBound:maxBound));
        
        %guess frequency
        [xmax,imax]=extrema(y);
        imax=sort(imax(1:5));
        avdif=0;
        N=0;
        for(i=1:4)
            if(imax(i+1)-imax(i)>50)
                avdif=avdif+(x(imax(i+1))-x(imax(i))); N=N+1;
            end
        end
        avdif=avdif/N;
        f=1/avdif;
        
        
        %guess phase
        if(y(1)<-.5*A)
            UnitCircle=1;
        elseif(y(1)<0)
            UnitCircle=2;
        elseif(y(1)<.5*A)
            UnitCircle=3;
        else
            UnitCircle=4;
        end
        
        GoingUp=(y(1)<y(3));
        
        
        switch UnitCircle
            case 1
                if(GoingUp)
                    phase = 15*pi/8;
                    %WAAAAAOW
            
        
        
        
     %x=1:.01/pi:2*pi;
     %y=.5*sin(3.*x+pi);
     optimset('TolX',1e-10);
     optimset('MaxIter',1e+10);
     [bestcoeffs,fval,exitflag]=fminsearch(@fun,[A f 5*pi/6],[],x,y)
     yfit=bestcoeffs(1)*sin(bestcoeffs(2)*2*pi.*x+bestcoeffs(3));
     %Now compare y with yfit
     figure;
     plot(x,y,x,yfit);
     hold on
     %plot(x,10*sin(2*pi*53.*x+0),'color','red');
     bestcoeffs



     function out=fun(coeff,X,Y)
         a = coeff(1);
         b = coeff(2);
         c = coeff(3);
         Y_fun = a .* sin(b*2*pi.*X+c);
         DIFF = Y-Y_fun; 
         SQ_DIFF = DIFF.^2;

         out = sum(SQ_DIFF);
     end

end