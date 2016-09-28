% calculate the absolute magnitude of torroidal flow, then plot average
% profile with error bars?
function absFlow(vel,impacts,time,rangeT,rangeP)
    ind1 = find(time>1,1);
    ind2 = find(time>2,1);
    
    [ave1,stdev1] = addNans(vel(ind1:ind2,rangeT));
    [ave2,stdev2] = addNans(vel(ind1:ind2,rangeP));
    flow = ave1-ave2;
    figure;plot(impacts(rangeT),flow);title('flow Profile');
    
    vel(:,rangeT) = vel(:,rangeT) - vel(:,rangeP);
    vel(:,rangeP) = vel(:,rangeP).* nan;
%     ind1 = find(time>1,1);
%     ind2 = find(time>2,1);
    %prof = mean(vel(ind1:ind2,rangeT),1); % average in the time direction
    size(vel)
    [ave,stdev] = addNans(vel(ind1:ind2,rangeT));
    [X,Y] = meshgrid(time(ind1:ind2),impacts(rangeT));
    figure; surf(X,Y,vel(ind1:ind2,rangeT)'); shading interp; view([ 0 90]); colorbar;
    figure; errorbar(impacts(rangeT),ave./2,stdev./2,stdev./2); 
    title('Average Toroidal Flow');
end
