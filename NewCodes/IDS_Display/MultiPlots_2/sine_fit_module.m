%% Calculate Sine Fit to Data for Multiplot_2
% Generates amplitude and phase initial guess from an FFT
% Fitting itself is done through sine_fit, which using nlnfit, which in
% turn uses Levenberg-Marquardt fitting.
% Also calculated percentace of FFT power spectrum at the injector
% frequency and its higher harmonics.
% Also caluclates the error in phase in a WEIRD WAY. 
function [guess,param,saveDat,SigDev,RMS,RMS_ideal] = sine_fit_module(in, doubleplot,dat,n,i,saveDat)
 

        %%%%%% Upper Array %%%%%%%%%%%
        % Calculate sine fit for uppe array
        for j = 1: 1+in(n).doubleplot
            size(dat(in(n).line).vel);
            if plotType==1
                signal = dat(in(n).line).vel(:,doubleplot(j,i));
            elseif plotType == 2
                signal = dat(in(n).line).temp(:,doubleplot(j,i));
            end
            if exist('savSig','var');size(savSig)
            end
            size(signal)
            %savSig(1,i,n,:)=signal;
            Fsamp = 1/(mean(diff(dat(1).time.*(in(n).timeScale.*1e-3))));
            offset = nanmean(signal);
            amp = max(signal)-offset;
            freq = 14500;
            phase = pi/2;


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

            guess(i,:,n,j) = [offset,max(signal)-offset,pi/2,freq];

            [param(i,(2:5)+5*(j-1),n),data(1:length(dat(1).time),i)] = sine_fit( ...
                dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
                guess(i,:,n,j),0);
            % Use FFT Result

            if param(i,3+5*(j-1),n)<0 % 180degree phase
                param(i,3+5*(j-1),n)=-param(i,3+5*(j-1),n);
                param(i,4+5*(j-1),n)=param(i,4+5*(j-1),n)+pi;
                disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 1, Impact: ' num2str(dat(1).impacts(i))]);
            end
            param(i,3+5*(j-1),n)= param_fft(i,3+5*(j-1),n);


            % Calculate fft reconstruction validity
            harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
            harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
            harm3=0;harm4=0;
            try harm3 = bandpower(signal-offset,Fsamp,[41,46].*1e3);end
            try harm4 = bandpower(signal-offset,Fsamp,[55.5,60.5].*1e3);end
            %harm5 = bandpower(signal-offset,Fsamp,[70,75].*1e3);
            ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal)) ]); % Nyquist
            pRel(i,j) = (harm1+harm2+harm3+harm4)/ptot;
            saveDat(n).FFT(i,j)=pRel(i,j);
            % calculate sine fit from FFT (again)

            % if the fit isnt valid, dont plot it
            try data((1:length(dat(1).time)+length(dat(1).time)*(j-1) ),...
                    i.*(pRel(i,j)<CutPow))=signal;end
            saveDat(n).FFT(i,j)=pRel(i,j);
            
            % attempt lm error analysis
             dp = [0.001, 0.001, 0.001, 0.001]; % fractional increment of 'p' for numerical derivatives
            [p_fit, Chi_sq, dPar(n,i,j,:), ~, corr, R2, cvg_hst] = ...
             lm(@SineFitLM, param(i,(2:5)+5*(j-1),n), dat(1).time.*(in(n).timeScale.*1e-3), signal', 0.001, dp);%, p_min,p_max,0)
            %pause(1);
%                         [BETA,R,J,COVB,MSE] = nlinfit(dat(1).time'.*(in(n).timeScale.*1e-3),signal',@SineFitNLN,param(i,2:5,n));
%                         parErrNLN(n,i,:)=sqrt(diag(sqrt(mean(R.^2))*pinv(J'*J)));
%                         SigDev(n,i,1)=parErrNLN(n,i,3);
%                         clear BETA R K COVB MSE
            %error('HALTING');
            % Calculate RMS Error
            % Calculate RMS Errorfor rmsPhase = 1:200
                 for rmsPhase = 1:200
                     %RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/10)))).^2));
                     RMS(j,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n) + (-pi +pi*rmsPhase/100))) ).^2 ));
                     RMS_ideal(j,i,n,rmsPhase)= sqrt(mean( ((param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n)))...
                         -(param(i,2+5*(j-1),n)+param(i,3+5*(j-1),n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4+5*(j-1),n) + (-pi +pi*rmsPhase/100)) ) ).^2 ));
                 end
%                             [p_fit_dat(n,i,2,:), Chi_sq, dPar_dat(n,i,2,:), ~, corr, R2, cvg_hst] = ...
%                             lm(@singletGauss1DLM, [-(RMS(2,i,n,1)-RMS(2,i,n,10))*pi,0,1.5,RMS(2,i,n,1)], (-pi +pi*(1:20)/10)', RMS(2,i,n,:), 0.0001, ones(1,4).*.001);%, p_min,p_max,0)
%                              p_fit_dat(n,i,1,:)=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(1,i,n,:))',@singletGauss1D,[-(RMS(1,i,n,1)-RMS(1,i,n,50))*pi,.5,1.5,RMS(1,i,n,1)]);
%                              [Y,I] = min( (squeeze(RMS_ideal(1,i,n,:))- RMS_ideal(1,i,n,1)/2).^2);% Find HWHM
%                              SigDev(n,i,1) = p_fit_dat(n,i,1,3)-abs(-pi +pi*(I)/100);
                 [p_fit_dat(n,i,j,:),R(n,i,j,:),J,COVB(n,i,j,:,:),MSE(n,i,j),ERRORMODELINFO(n,i,j)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(j,i,n,:))',@singletGauss1D,[-(RMS(j,i,n,1)-RMS(1,i,n,100))*pi,.5,1.5,RMS(j,i,n,1)]);
                 [p_fit_dat_ideal(n,i,j,:),R_ideal(n,i,j,:),J,COVB_ideal(n,i,j,:,:),MSE_ideal(n,i,j),ERRORMODELINFO_ideal(n,i,j)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS_ideal(j,i,n,:))',@singletGauss1D,[-(RMS_ideal(j,i,n,1)-RMS_ideal(1,i,n,100))*pi,.5,1.5,RMS_ideal(j,i,n,1)]);

%                               [Y,I] = min( (squeeze(RMS_ideal(1,i,n,:))- RMS_ideal(1,i,n,1)/2).^2);% Find HWHM
%                               SigDev(n,i,1) = p_fit_dat(n,i,1,3)-abs(-pi +pi*(I)/100);
                    SigDev(n,i,j) = abs(p_fit_dat(n,i,j,3)-p_fit_dat_ideal(n,i,j,3)); % Delta sigmas is error
        end
%             %%%%%% Lower Array %%%%%%%%%%%%%%%
%             if in(n).doubleplot
%                 % Fit Sine
%                 if plotType==1
%                     signal = dat(in(n).line).vel(:,doubleplot(2,i));
%                 elseif plotType == 2
%                     signal = dat(in(n).line).temp(:,doubleplot(2,i));
%                 end
%                 %savSig(2,i,n,:)=signal;
%                 offset = mean(signal);
%                 amp = max(signal)-offset;
% 
%                 xfft(:,i,n,2)=fft(signal-offset);
%                 P1=xfft(1:floor(length(signal)/2)+1,i,n,2);
%                 P2=abs(P1);
%                 P2(2:end-1)=2*P2(2:end-1);
%                 f=Fsamp*(0:length(signal)-1)/length(signal);
%                 [Y,I]=max(P2);
%                 amp = P2(I)/length(signal);
%                 phase = mod(unwrap(angle(P1(I))),2*pi);
%                 data_fft(1:length(dat(1).time),i,n)=offset+amp*sin(f(I)*2*pi*dat(1).time.*(in(n).timeScale.*1e-3) + phase);
%                 param_fft(i,7:10,n) = [offset, amp, phase, f(I)];
%                 if f(I)>20000 ; f(I)=f(I)/2; end % hit harmoinc
% 
%                 guess(i,:,n,2) = [offset,max(signal)-offset,pi/2,freq];
%                 [param(i,7:10,n),data(length(dat(1).time)+1:2*length(dat(1).time),i)] = ...
%                     sine_fit(dat(1).time'.*(in(n).timeScale.*1e-3),signal',[nan,nan,nan,freq], ...
%                     guess(i,:,n,2),0);
%                 %pause(1);
%                 % Use the FFT Results
%                 %param(i,8,n)= abs(param_fft(i,8,n));
%                 if param(i,8,n)<0 % 180degree phase
%                     param(i,8,n)=-param(i,8,n);
%                     param(i,9,n)=param(i,9,n)+pi;
%                     disp([' WARNING: NEGATIVE AMPLITUDE @ n=' num2str(n) ', Line 2, Impact: ' num2str(dat(1).impacts(i))]);
%                 end
%                 param(i,8,n)= (param_fft(i,8,n)); % needs to go here, otherwise sine of amp cant get fixed
%                 % Calculate fft reconstruction validity
%                 harm1 = bandpower(signal-offset,Fsamp,[13,17].*1e3);
%                 harm2 = bandpower(signal-offset,Fsamp,[26,31].*1e3);
%                 harm3=0;harm4=0;
%                 try harm3 = bandpower(signal-offset,Fsamp,[41,46].*1e3);end
%                 try harm4 = bandpower(signal-offset,Fsamp,[55.5,60.5].*1e3);end
%                 ptot =  bandpower(signal-offset,Fsamp,[0,Fsamp*(length(signal)-2)/(2*length(signal))]);
%                 pRel(i,2) = (harm1+harm2+harm3+harm4)/ptot;
%                 xfft(:,i,n,2)=abs(fft(signal-offset))/length(signal);
%                 % if the fit isnt valid, dont plot it
%                 try data(length(dat(1).time)+1:2*length(dat(1).time)...
%                         ,i.*(pRel(i,2)<CutPow))=signal;end
%                 saveDat(n).FFT(i,2)=pRel(i,2);
% 
%                 % attempt lm error analysis
% %                             dp = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]; % fractional increment of 'p' for numerical derivatives
% %                             [p_fit, Chi_sq, dPar(n,i,1,:), ~, corr, R2, cvg_hst] = ...
% %                             lm(@SineFitLM,  param(i,7:10,n), dat(1).time.*(in(n).timeScale.*1e-3), signal, 0.0001, dp);%, p_min,p_max,0)
% %                              [BETA,R,J,COVB,MSE] = nlinfit(dat(1).time'.*(in(n).timeScale.*1e-3),signal',@SineFitNLN,param(i,7:10,n));
% %                             parErrNLN(n,i,:)=sqrt(diag(sqrt(mean(R.^2))*pinv(J'*J)));
% %                             SigDev(n,i,2)=parErrNLN(n,i,3);
% %                             clear BETA R K COVB MSE
% 
%                 % Calculate RMS Errorfor rmsPhase = 1:200
%                  for rmsPhase = 1:200
%                      %RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,2,n)+param(i,3,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,4,n) + (-pi +pi*rmsPhase/10)))).^2));
%                      RMS(2,i,n,rmsPhase) = sqrt(mean( (signal'-(param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n) + (-pi +pi*rmsPhase/100))) ).^2 ));
%                      RMS_ideal(2,i,n,rmsPhase)= sqrt(mean( ((param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n)))-(param(i,7,n)+param(i,8,n)*sin(2*pi*dat(1).time'.*(in(n).timeScale).*1e-3.*14500 + param(i,9,n) + (-pi +pi*rmsPhase/100)) ) ).^2 ));
%                  end
% %                             [p_fit_dat(n,i,2,:), Chi_sq, dPar_dat(n,i,2,:), ~, corr, R2, cvg_hst] = ...
% %                             lm(@singletGauss1DLM, [-(RMS(2,i,n,1)-RMS(2,i,n,10))*pi,0,1.5,RMS(2,i,n,1)], (-pi +pi*(1:20)/10)', RMS(2,i,n,:), 0.0001, ones(1,4).*.001);%, p_min,p_max,0)
%                  [p_fit_dat(n,i,2,:),R(n,i,2,:),J,COVB(n,i,2,:,:),MSE(n,i,2),ERRORMODELINFO(n,i,2)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS(2,i,n,:))',@singletGauss1D,[-(RMS(2,i,n,1)-RMS(2,i,n,100))*pi,.5,1.5,RMS(2,i,n,1)]);
%                  [p_fit_dat_ideal(n,i,2,:),R_ideal(n,i,2,:),J,COVB_ideal(n,i,2,:,:),MSE_ideal(n,i,2),ERRORMODELINFO_ideal(n,i,2)]=nlinfit((-pi +pi*(1:200)/100),squeeze(RMS_ideal(2,i,n,:))',@singletGauss1D,[-(RMS_ideal(2,i,n,1)-RMS_ideal(2,i,n,100))*pi,.5,1.5,RMS_ideal(2,i,n,1)]);
% 
% %                              [Y,I] = min( (squeeze(RMS_ideal(2,i,n,:))- RMS_ideal(2,i,n,1)/2).^2);% Find HWHM
% %                              SigDev(n,i,2) = p_fit_dat(n,i,2,3)-abs(-pi +pi*(I)/100);
%                     SigDev(n,i,2) = abs(p_fit_dat(n,i,2,3)-p_fit_dat_ideal(n,i,2,3)); % Delta sigmas is error
%             end
% 
% 
%         %catch
%         %    display(['halted at ' num2str(i)]);
%         %end
end


