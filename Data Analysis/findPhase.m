% find all the phases

data = dat(1).vel;
time = dat(1).time;

param = zeros(15,5);

for i = 1:15
    k = 131;
    while ~isnan(data(k,7+i)) %&& k>69
        k = k -1;
    end
    param(i,5) = k+1;
    param(i,1:4) = lsqcurvefit(@fitSine,[.3726,20515,.1801,.5209],time((k+1):131)*1e-3,abs(diff(data((k+1):132,(7+i)))),[0,10000,0,0],[3,40000,2,1]);
    figure; plot(time((k+1):131),abs(diff(data((k+1):132,7+i))),time((k+1):131),fitSine(param(i,1:4),time((k+1):131)*1e-3));
end


param1 = lsqcurvefit(@fitSine,[2,20000,1,0],dat(1).ItorTime(1732:2028)*1e-3,dat(1).Itor(1732:2028)-mean(dat(1).Itor(1732:2028)),[0,10000,0,0],[3,40000,2,1]);
figure; plot(dat(1).ItorTime(1732:2028),dat(1).Itor(1732:2028)-mean(dat(1).Itor(1732:2028)), dat(1).ItorTime(1732:2028),fitSine(param1,dat(1).ItorTime(1732:2028)*1e-3));