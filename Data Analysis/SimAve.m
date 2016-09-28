% frame averaging to match NIDS

NIDS = importdata('dat815031802810.mat');
IDS = importdata('dat15040102710.mat');

Pts = 7:11;
IDS.time = IDS.time.*1e-6;

slidefactor = 1.07;
NIDS.time = NIDS.time + slidefactor;

%remove NaNs
for passes =  1:3 % go through three times to make sure
    for m = 12:28 % chans
        for n = 2:(length(NIDS.time)-1) % timepts
            if isnan(NIDS.vel(n,m))
                assignin('base','n',n);
                assignin('base','m',m);
                NIDS.vel(n,m)=mean([NIDS.vel(n-1,m),NIDS.vel(n+1,m)]);
            end
        end
    end
end


aveVel = zeros(Pts(end),28);
x=1;

for i = Pts
    counter = 1;
    while(x<= 105 && NIDS.time(x)<IDS.time(i))
        aveVel(i,:) = aveVel(i,:)+NIDS.vel(x,:);
        counter = counter +1;
        x=x+1;
    end
    aveVel(i,:) = aveVel(i,:)./counter;
end
    