%  find PIX_SP by fitting to known spectral lines
close all;
% data
shot = 150625140;


% wavelength
lambda = 1e-9 * [464.913, 464.742, 464.308, 464.181]; % descending order


trimTime = [1];
trimChan = [61:343];
addpath('S:\MC_IDS\Matlab Code\Core Codes\lsqcurvefit');
addpath('T:\PhantomMovies');
data = importdata(['shot' int2str(shot) '.mat']);
data = squeeze((sum(data(trimTime,:,trimChan),1)));
data = squeeze(sum(data,2));
data = data(end:-1:1);
% find peaks
%  addpath('T:\RChandra\A-A-Ron Code\General Matlab\extrema');
%  [~,imax,~,~] = extrema(data);
% find the peaks by cheating
peaks = [78,63,26,15];
plot(data,'k')
hold on;

for i = 1:4
    guess(i,:) = [ data(peaks(i))*1.5,peaks(i), 2, min(data( (peaks(i)-7):(peaks(i)+7) ))];
    PreParam(i,:) = lsqcurvefit(@singletFun, guess(i,:), [(peaks(i)-7):(peaks(i)+7) ]', data( (peaks(i)-7):(peaks(i)+7) ));
    plot([(peaks(i)-7):(peaks(i)+7) ]', singletFun(PreParam(i,:),[(peaks(i)-7):(peaks(i)+7) ]'),'-*b');
    plot([(peaks(i)-7):(peaks(i)+7) ]', singletFun(guess(i,:),[(peaks(i)-7):(peaks(i)+7) ]'),'-*r');
end

figure;plot(PreParam(end:-1:1,2)',lambda(end:-1:1),'*r');
p = polyfit(PreParam(end:-1:1,2)',lambda(end:-1:1),1);
p(1)*1e6
hold on;
plot(PreParam(end:-1:1,2)',[p(2)+p(1)*PreParam(end:-1:1,2)']);