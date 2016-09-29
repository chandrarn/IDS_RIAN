
% Add lots of shots together
%shots = [ 150625173:150625203,150625205:150625223];%,150625257:150625283];
%shots = [150625257:150625283];
%shots = [150626126:150626137]
% shots = [ 150624149:150624163];
% shots = [150626174:150626193];
%shots = [ 150707010:150707017];
%shots = [151022102:151022109,151022111:151022129];
%shots = 151027102:151027136
%shots = [151022115,151022108,151022114,151022117,151022111,151022121];
shots = 151022108;
shotnum = [151022020]; % save shot number

cd('T:\PhantomMovies');
data = zeros(484,96,352);
t = zeros(484,1);
% MetaShot = zeros(length(shots), 30,96,352);
% cd('H:\Cines');
% shots = [ 150630115:150630315];
% data = zeros(96,352,527);
% t = zeros(527,1);
ind = 1;
for i = shots
    try
        temp = double(importdata(['shot' num2str(i) '.mat']));%importdata(['Shot ' num2str(i) '.mat']);%
%         for a = 1:size(data,1)
%             for j = 1:size(data,2)
%                 for k = 1:size(data,3)
%                     data(a,j,k) = data(a,j,k) + temp(a,j,k);%temp.CineArray;
%                 end
%             end
%          end
        data = data + temp;
        t = t +importdata(['t' num2str(i) '.mat']);% temp.TimeVector; %
    %     MetaShot(ind,:,:,:) = temp;
    %     ind = ind +1
        i
    catch
        display(['IMPORT FAILED, SHOT: ' num2str(i)]);
    end
end
data = data ./length(shots);
t = t./length(shots);

if ~isempty(shotnum)
     save(['shot' num2str(shotnum) '.mat'],'data');
     save(['t' num2str(shotnum) '.mat'],'t');
end
% stdBlock = zeros(30,96,352);
% for i = 1:size(data,1)
%     for j = 1:size(data,2)
%         for k = 1:size(data,3)
%             stdBlock(i,j,k) = std(MetaShot(:,i,j,k));
%         end
%     end
% end

