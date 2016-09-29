function x = sihi_smooth(y,time,inj_freq)
% take a boxcar average over an injector cycle
% interpolate data to uniform time base

    %Start with check in case using imaginary numbers:
    doImag = 0;
    for iy = 1:length(y)
        if imag(y(iy)) ~= 0
            doImag = 1;
            break
        end
    end
    if doImag == 1
        yi = imag(y);
        yr = real(y);
        injCyc = 1/inj_freq;
        Navg = 100;
        tint = time(1):injCyc/Navg:time(end); %make it 100 time points per injector cycle
        yrint = interp1(time,yr,tint);
        yiint = interp1(time,yi,tint);
        %
        xrint = zeros(1,length((Navg/2+1):length(tint)-(Navg/2)));
        xiint = zeros(1,length((Navg/2+1):length(tint)-(Navg/2)));
        for it = (Navg/2+1):length(tint)-(Navg/2)
            xrint(it) = mean(yrint(it-(Navg/2):it+Navg/2-1));
            xiint(it) = mean(yiint(it-(Navg/2):it+Navg/2-1));
        end
        xiint(1:Navg/2) = xiint(Navg/2+1);
        xiint(length(tint)-Navg/2+1:length(tint)) = xiint(length(tint)-Navg/2);
        xrint(1:Navg/2) = xrint(Navg/2+1);
        xrint(length(tint)-Navg/2+1:length(tint)) = xrint(length(tint)-Navg/2);
        xr = interp1(tint,xrint,time);
        xi = interp1(tint,xiint,time);
        x = xr + sqrt(-1)*xi;
    else
        injCyc = 1/inj_freq;
        Navg = 100;
        tint = time(1):injCyc/Navg:time(end); %make it 100 time points per injector cycle
        yint = interp1(time,y,tint);
        %
        xint = zeros(1,length((Navg/2+1):length(tint)-(Navg/2)));
        for it = (Navg/2+1):length(tint)-(Navg/2)
            xint(it) = mean(yint(it-(Navg/2):it+Navg/2-1));
        end
        xint(1:Navg/2) = xint(Navg/2+1);
        xint(length(tint)-Navg/2+1:length(tint)) = xint(length(tint)-Navg/2);
        x = interp1(tint,xint,time);
    end
    x(end) = x(end-1);
return
end