% Display all the current data for some given shot
% Rian Chandra

Shots = [129810, 129213, 128585];

cd('T:\IDS\Data Repository');
for i= 1:length(Shots);
    shot(i)=eval(sprintf('load(''dat%0.0f'');', i));
end


