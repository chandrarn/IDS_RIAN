function [dbrdj, dbpdj, dbtdj] = ...
    imp_pol_deriv(top, mid, bot, array, order)

% this takes the poloidal derivative using data from the IMP
% using t, m, b specify the main stem of the calculation
% 1t is top - mid
% 1b is bot - mid
% 2 is the average of 1t and 1b

% b is the signal of interest
% N is the number of probes

% the tor and pol distance between the mid and either 
% top or bot array is 0.875 inches
dr = 0.0222;

if array == 't'
    dbdj = (top - mid)/dr;
    
elseif array == 'm'
    if(strcmp(order, '1t'))
        dbdj = (top - mid)/dr;
    elseif(strcmp(order, '1b'))
        dbdj = (bot - mid)/dr;
    elseif(strcmp(order, '2'))
        dbdj = (top - 2*mid + bot)/(2*dr);
    end
    
elseif array == 'b'
    dbdj = (bot - mid)/dr;
    
end

dbrdj = squeeze(dbdj(1, :, :));
dbpdj = squeeze(dbdj(2, :, :));
dbtdj = squeeze(dbdj(3, :, :));


