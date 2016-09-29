function[n_start, n_end, n_time] = shotDuration(shot, time, t_override)
% Function whose purpose is to determine the duration of the plasma

if t_override == 0
    decay_t = 0.6; % time to continue analyzing after SPAs have turned off
                      % (ie: current decay time [ms])

    spa_length = 6; % Length of SPA_ENABLE as presently set up.  If the length
                    % is different from this, the form of this array has
                    % changed and this code will throw an error.

    mdsopen('landau.hit::hitsi', shot);
    spa_x = mdsvalue('dim_of(\SPA_ENABLE_X)');
    spa_y = mdsvalue('dim_of(\SPA_ENABLE_Y)');
    mdsclose();

    if or(length(spa_x) ~= spa_length, length(spa_y) ~= spa_length)
        disp('The form of SPA_ENABLE_X,Y has changed and may cause this function to set the shot length incorrectly');
        n_start = 1; % use entire movie
        n_end = length(time);
        n_time = length(time);
        return;
    end

    t_start = 0;
    t_end = decay_t + max(1e3*[spa_x(5) spa_y(5)]);
else
    t_start = t_override(1);
    t_end = t_override(2);
end

n_start = find(time >= t_start, 1);
n_end = find(time >= t_end, 1);
n_time = n_end - n_start + 1;

end