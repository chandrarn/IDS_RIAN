
% plot spectrum with wavelengths
%shot = 151217026;
shot = 160518012;
%try
%    temp = CineArray(1,1,end:-1:1);
%    temp = dat(1).raw(1,1,end:-1:1);
%catch
    addpath('T:\PhantomMovies');
    load(['Shot ' num2str(shot) '.mat']);
    load(['dat' num2str(shot) '.mat']);
    %CineArray = CineArray(end:-1:1,:,:);
%end

Nist = { 464.181, 'OII'; 464.308, 'NII'; 464.313, 'WI';  464.339 'OII'; 464.616, 'CrI'; 464.681, 'FeII';...
    464.718, 'FeI'; 464.742, 'CIII'; 464.78, 'OII';...
    464.861, 'AlII'; 464.887, 'CrI'; 464.913, 'OII'; 465.025, 'CIII'; 465.084, 'OII'; 465.147, 'CIII'; ...
    465.3,'AlII'; };


%spectrum = squeeze(sum(sum(CineArray,2),3));
spectrum = sum(CineArray(end:-1:1,:,:),2);
figure('name','Spectrum Plot','position',[5,580,1670,400;]);
wvlngt = ((1:length(spectrum))-mean(dat(1).param.Center(:,1)) -3-.5)*1.153e-11 +dat(1).param.LineLam(3); % replace 1 with 2

plot(wvlngt*1e9,spectrum./1e5,'-*','LineWidth',2);hold on;
xlabel('\lambda [nm]'); ylabel('Intensity [arb]');title('IDS Spectrum');
box on; grid on; 
%set(gca,'xtick',[464.4:.05:465.6]);
set(gca,'xtick',[464.0:.075:465.5]);
% Plot NIST
for i = 1:size(Nist,1)
    plot([Nist{i,1},Nist{i,1}],[2.8,3.8],'r');
    text(Nist{i,1}+.005,3.825,Nist{i,2},'fontsize',12,'fontweight','bold');
end