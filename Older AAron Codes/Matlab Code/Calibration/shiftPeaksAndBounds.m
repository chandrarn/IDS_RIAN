% Code to take original 'peaks' and 'bounds' data from calibration in the 
% tree and transform it to match new (presumably smaller) CCD size used for
% plasma data.

saveToTree = 1; % safety feature for testing
peaksRefShot = 126684; % reference shot from which to pull 'peaks', 'bounds', and 'center'
saveToShots = [126688:126693]; % shot list of plasma data with different CCD size

x_shift = 126; % shift in the horizontal direction
y_shift = 129; % shift in the vertical direction

% 126, 129 for 512x384 -> 256x128

mdsopen('landau.hit::analysis', peaksRefShot);
peaks = mdsvalue('\IDS_PEAKS');
bounds = mdsvalue('\IDS_BOUNDS');
center = mdsvalue('\IDS_CENTER');
mdsclose();

peaks(:, 2) = peaks(:, 2) - x_shift;
peaks(:, 3) = peaks(:, 3) - y_shift;
bounds = bounds - y_shift;
center = center - y_shift;

if saveToTree
    for n = 1:length(saveToShots)
        mdsopen('landau.hit::analysis', saveToShots(n));
        
        disp(['Editing shot: ' num2str(saveToShots(n))]);
        mdsput('\IDS_PEAKS', '$', peaks)
        mdsput('\IDS_BOUNDS', '$', bounds)
        mdsput('\IDS_CENTER', '$', center)
        
        mdsclose();
    end
else
    saveToShots'
    peaks
    bounds
    center
end