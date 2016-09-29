%addnans
function [output,stdev] = addNans(data) %updated to average anything it's given

%data = dat.vel;
output = zeros(1,size(data,2));
temp = zeros(1,size(data,1));
stdev = zeros(1,size(data,2));
for i = 1:size(data,2)
    counter = 1;
    for k = 1:size(data,1)
        if ~isnan(data(k,i))
            output(i) = output(i)+data(k,i);
            temp(counter) = data(k,i);
            counter = counter +1;
        end
    end
    stdev(i) = std(temp(1:counter-1));
    output(i) = output(i)/counter;
end

%figure;plot(output);
        
data1 = data;
for i = 1:size(data,2)
    data1(:,i) = data1(:,i)-output(i);
end
%figure; surf(data1); shading interp; view([ 0 90]);colorbar

assignin('base','data1',data1);
end