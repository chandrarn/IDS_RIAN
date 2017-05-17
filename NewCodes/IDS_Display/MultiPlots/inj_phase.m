%% Find Injector Phase for Multiplot_2

function [injParam] = inj_phase(dat,timebound,in)

if n==1
    clear signal
    freq = 14500;
    
    if (in(1).shot <129500)|| (in(1).shot == 8129499)
        [Y,I1]=min( (dat(1).iinjxTime - timebound(1) ).^2 );
        [Y,I2]=min( (dat(1).iinjxTime - timebound(2) ).^2 );
        signal(:,1) = dat(1).iinjx(I1:I2);
        signal(:,2) = dat(1).iinjy(I1:I2);
        for i = 1:2
            offset = mean(signal(:,i));
            amp = max(signal(:,i))-offset;
            [injParam(i,:),~] = ...
                sine_fit(double(dat(1).iinjxTime(I1:I2))'.*1e-3,double(signal(:,i))',[nan,nan,nan,freq], ...
                [offset,amp,0,freq],0);
        end
    else
        [Y,I1]=min( (dat(1).iinjaTime - timebound(1) ).^2 );
        [Y,I2]=min( (dat(1).iinjaTime - timebound(2) ).^2 );
        signal(:,1) = dat(1).iinja(I1:I2);
        signal(:,2) = dat(1).iinjb(I1:I2);
        signal(:,3) = dat(1).iinjc(I1:I2);
        offset = mean(signal);
        amp = max(signal)-offset;
        for i = 1:3
            offset = mean(signal(:,i));
            amp = max(signal(:,i))-offset;
            [injParam(i,:),~] = ...
                sine_fit(double(dat(1).iinjaTime(I1:I2))'.*1e-3,double(signal(:,i))',[nan,nan,nan,freq], ...
                [offset,amp,0,freq],0);
        end
    end
    for i = 1:size(injParam,1)
        if injParam(i,2)<0 % 180degree phase
            injParam(i,2)=-injParam(i,2);
            injParam(i,3)=injParam(i,3)+pi;
            disp(['WARNING: INJECTOR ' num2str(i) ' NEGATIVE AMPLITUDE']);
        end
    end

end

end