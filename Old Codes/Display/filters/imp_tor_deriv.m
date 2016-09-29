function [dbrdk, dbpdk, dbtdk] = ...
    imp_tor_deriv(top, mid, bot, array, order)

% this takes the poloidal derivative using data from the IMP

% b is the signal of interest
% N is the number of probes

% the tor and pol distance between the mid and either 
% top or bot array is 0.875 inches
dr = 0.0222;

if array == 't'
    if(strcmp(order, '1'))
        dbdk = (mid - top)/dr;
    elseif(strcmp(order, '2'))
        dbdk = (bot - top)/(2*dr);
    end

elseif array == 'm'
    if(strcmp(order, '1t'))
        dbdk = (mid - top)/dr;
    elseif(strcmp(order, '1b'))
        dbdk = (bot - mid)/dr;
    elseif(strcmp(order, '2'))
        dbdk = (bot - top)/(2*dr);
    end

elseif array == 'b'
    if(strcmp(order, '1'))
        dbdk = (bot - mid)/dr;
    elseif(strcmp(order, '2'))
        dbdk = (bot - top)/(2*dr);
    end
    
end

dbrdk = squeeze(dbdk(1, :, :));
dbpdk = squeeze(dbdk(2, :, :));
dbtdk = squeeze(dbdk(3, :, :));

