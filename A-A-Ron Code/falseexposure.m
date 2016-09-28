data1 = data(1:210,:,:);

temp = zeros(21,96,320);
for i = 0:20
    for j = 1:10
        for l = 1:96
            for m = 1:320
                temp(i+1,l,m) = temp(i+1,l,m) + data1(i*10 + j,l,m);
            end
        end
    end
end

temp = temp./10;
