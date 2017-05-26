%% Calculate Sine Fit to Data for Multiplot_2
% Generates amplitude and phase initial guess from an FFT
% Fitting itself is done through sine_fit, which using nlnfit, which in
% turn uses Levenberg-Marquardt fitting.
% Also calculated percentace of FFT power spectrum at the injector
% frequency and its higher harmonics.
% Also caluclates the error in phase in a WEIRD WAY. 

% PARAM: [IMPACT, OFFSET, AMP, PHASE, FREQ]
function [guess_out,param_out,saveDat,SigDev_out,RMS_out,RMS_ideal_out,data_out,pRel_out,dPar_out] = ...
    sine_fit_module(in, doubleplot,dat,n,i,saveDat,plt)

        for j = 1: 1+in(n).doubleplot % loop through arrays to calc fit
            % Extract signal to fit to
            if plt.Type==1
                signal = dat(in(n).line).vel(:,doubleplot(j,i));
            elseif plotType == 2
                signal = dat(in(n).line).temp(:,doubleplot(j,i));
            end
            
            % Set up initial guesses
            Fsamp = 1/(mean(diff(dat(1).time.*(in(n).timeScale.*1e-3))));
            offset = nanmean(signal);
            amp = max(signal)-offset;
            freq = 14500;
            phase = pi/2;

            % Perform FFT
            xfft(:,i,n,j)=fft(signal-offset);
            P1=xfft(1:floor(length(signal)/2)+1,i,n,j);
            P2=abs(P1);
            P2(2:end-1)=2*P2(2:end-1);
            f=Fsamp*(0:length(signal)-1)/length(signal);
            [Y,I]=max(P2);
            amp = P2(I)/length(signal);
            phase = mod(unwrap(angle(P1(I))),2*pi);
            data_fft(1:length(dat(1).time),i,n)=offset+amp*sin(f(I)*2*pi*dat(1).time.*(in(n).timeScale.*1e-3) + phase);
            param_fft(i,(2:5)+5*(j-1),n) = [offset, amp, phase, f(I)];
            if f(I)>20000 ; f(I)=f(I)/2; end % hit harmoinc

            guess(i,:,n,j) = [offset,max(signal)-offset,pi/2,freq]; % save guesses
            
            % Perform LM sine function fit, with FFT guesses
            [param(i,(2:5)+5*(j-1),n),data((1:length(dat(1).time))...
                +length(dat(1).time)*(j-1),i)] = sine_fit( ...
                dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                guess(i,:,n,j),0);
            
            % Sometimes the amplitude of the fit is negative. +Pi phase
            % shift to compinsate
            if param(i,3+5*(j-1),n)<0 % 180degree phase
                param(i,3+5*(j-1),n)=-param(i,3+5*(j-1),n);
                param(i,4+5*(j-1),n)=param(i,4+5*(j-1),n)+pi;
                disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 1, Impact: ' num2str(dat(1).impacts(i))]);
            end
            
            param(i,3+5*(j-1),n)= param_fft(i,3+5*(j-1),n);


            % Calculate fft reconstruction validity: What fraction of the
            % total spectral power is the injector frequency? 
            % The window around the injector frequency and its harmonics
            % accounts for FFT sampling: it may not find the injectors at
            % EXACTLY the Inj frequency, power may be spread over a few
            % frequency bins.
            harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
            harm2=0;harm3=0;harm4=0;
            try harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);end 
            try harm3 = bandpower(signal-offset,Fsamp,[41,46].*1e3);end
            try harm4 = bandpower(signal-offset,Fsamp,[55.5,60.5].*1e3);end
            %harm5 = bandpower(signal-offset,Fsamp,[70,75].*1e3);
            ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal)) ]); % Nyquist
            pRel(i,j) = (harm1+harm2+harm3+harm4)/ptot;
            saveDat(n).FFT(i,j)=pRel(i,j);
            

            % if the fit isnt valid, dont plot it
            try data((1:length(dat(1).time)+length(dat(1).time)*(j-1) ),...
                    i.*(pRel(i,j)<CutPow))=signal;end
            saveDat(n).FFT(i,j)=pRel(i,j);
            
            % attempt lm error analysis
             dp = [0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives
            [p_fit, Chi_sq, dPar(n,i,j,:), ~, corr, R2, cvg_hst] = ...
             lm(@SineFitLM, param(i,(2:5)+5*(j-1),n), dat(1).time.*(in(n).timeScale.*1e-3), signal', 0.001, dp);%, p_min,p_max,0)
           
            % Calculate RMS Error
            % Calculate RMS Errorfor rmsPhase = 1:200
                 for rmsPhase = 1:200
                     %RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/10)))).^2));
                     RMS(j,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n) + (-pi +pi*rmsPhase/100))) ).^2 ));
                     RMS_ideal(j,i,n,rmsPhase)= sqrt(mean( ((param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n)))...
                         -(param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n) + (-pi +pi*rmsPhase/100)) ) ).^2 ));
                 end                           
                 [p_fit_dat(n,i,j,:),R(n,i,j,:),J,COVB(n,i,j,:,:),MSE(n,i,j),ERRORMODELINFO(n,i,j)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(j,i,n,:))',@singletGauss1D,[-(RMS(j,i,n,1)-RMS(1,i,n,100))*pi,.5,1.5,RMS(j,i,n,1)]);
                 [p_fit_dat_ideal(n,i,j,:),R_ideal(n,i,j,:),J,COVB_ideal(n,i,j,:,:),MSE_ideal(n,i,j),ERRORMODELINFO_ideal(n,i,j)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS_ideal(j,i,n,:))',@singletGauss1D,[-(RMS_ideal(j,i,n,1)-RMS_ideal(1,i,n,100))*pi,.5,1.5,RMS_ideal(j,i,n,1)]);

                    SigDev(n,i,j) = abs(p_fit_dat(n,i,j,3)-p_fit_dat_ideal(n,i,j,3)); % Delta sigmas is error
        end

        % Output Stuff
        guess_out=guess(i,:,n,:);
        param_out=param(i,:,n);
        SigDev_out=SigDev(n,i,:);
        RMS_out=RMS(:,i,n,:);
        RMS_ideal_out=RMS_ideal(:,i,n,:);
        data_out=data(:,i);
        pRel_out=pRel(i,:);
        dPar_out=dPar(n,i,:,:);

end


