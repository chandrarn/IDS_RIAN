% find phase for PDC shots
import MDSplus.*
Conn = Connection('landau.hit');

shots = [151022102:151022109, 151022111:151022128];

ref = 15102202210;
try
    dat;
catch
    load(['T:\IDS\Data Repository\dat' num2str(ref) '.mat']);
end
close all
figure;
Conn.openTree('pdc3',151022102);
Iinja = NATIVEvalue(Conn.get('\I_INJ_A'));
IinjT = double(NATIVEvalue(Conn.get('\I_INJ_A').getDimensionAt(0))).*1e-4;
%plot(dat(1).iinjaTime(17000:17200),dat(1).iinja(17000:17200)); hold on;
plot(IinjT(17000:17200),Iinja(17000:17200).*1e-3); hold on;
[~,n1] = min(abs(IinjT -1.7732)); %reference index closest to first peak
[~,n2] = min(abs(IinjT -1.7822)); %reference index closest to first trough
[~,n3] = min(abs(IinjT -1.7910)); %reference index closest to second peak
%plot([dat(1).iinjT(n2),dat(1).iinja(n2)],[-5000,5000],'g*');
n2 = 184;
for i = shots
        load(['T:\PhantomMovies\t' num2str(i) '.mat']);
        %Conn.openTree('pdc3',i);
        %iinja = NATIVEvalue(Conn.get('\I_INJ_A'));
        %iinjT = double(NATIVEvalue(Conn.get('\I_INJ_A').getDimensionAt(0))).*1e-4;
        
        %pick the amplitude at index n2, arbritrailly cooresponding to
        %trough. Map new point into graph relative to this point
%         if iinja(n2) > iinja(n2+2) % we're going down
%             lim1 = n1;
%             lim2 = n2;
%         else % going up
%             lim1 = n2;
%             lim2 = n3; 
%         end
%         
        %[~,Tloc] = min(abs(Iinja(lim1:lim2) - iinja(n2)));
%         figure;plot(iinjT(17000:17200),iinja(17000:17200))
%         Tloc = Tloc +lim1;
        %plot([IinjT(Tloc),IinjT(Tloc)],[-5,5],'r-*');
        plot([t(n2),t(n2)],[-5,5],'r-*');
        plot([t(n2+1),t(n2+1)],[-5,5],'g-*');
        X(find(shots == i),1) = t(n2);
        X(find(shots == i),2) = i;
        display(['Shot: ' num2str(i) ', t = ' num2str(t(n2))]);
end

[T(:,1),I] = sort(X(:,1));
T(:,2) = X(I,2);