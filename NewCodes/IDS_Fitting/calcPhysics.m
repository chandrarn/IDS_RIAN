function[temp, vel, int, param] = calcPhysics(fit_par, param, s, nn)

[n_time, n_chan, n_par] = size(fit_par); % determine size of data array

%% Temperature
% Calculate Instrument Temperatures

sigym = param.PIX_SP .* param.peaks(:, 5);
param.Inst_Temp(:, nn) = sigym.^2 * param.c^2 * param.IonMass(nn) / ...
    (param.CalLam^2 * param.kBoltz * 11605);

% Convert Instrument Temperature to dimensions of data

Inst_Temp = meshgrid(param.Inst_Temp(:, nn), 1:n_time);

% Expand 'PIX_SP' into same dimensions as Guassian width data

PIX_SP = meshgrid(param.PIX_SP, 1:n_time);

% Calculate Temperature

UnCorr_Temp = squeeze(fit_par(:, :, 5)).^2 .* PIX_SP.^2 .* param.c^2 .* param.IonMass(nn) ./ ...
    (mean(param.LineLam(nn))^2 * param.kBoltz * 11605);

temp = UnCorr_Temp - Inst_Temp; % subtract instrument temperature

%% Velocity
% Expand 'Center' into same dimensions as Gaussian offset data

center = meshgrid(param.Center(:, nn), 1:n_time);

% calculate velocity

vel = param.c * (squeeze(fit_par(:, :, 3)) - center) .* PIX_SP ./ param.LineLam(nn);

%% Intensity

% convert volume to amplitude

amp = squeeze(fit_par(:, :, 1)) ./ (2*pi * squeeze(fit_par(:, :, 4)) .* squeeze(fit_par(:, :, 5)));

% convert amplitude to area

area = sqrt(2*pi) * squeeze(fit_par(:, :, 5)) .* amp;

if isempty(s.sim) % real IDS data

    [rel_int, dummy] = meshgrid(param.REL_INT, 1:n_time); % expand REL_INT onto data sized grid
    size(rel_int)
    size(area)
    int = rel_int .* area;
    disp('Intensity corrected because data is real data from IDS');

else % Simulation Data - DO NOT correct for fiber throughput
    int = area;
    disp('Intensity NOT corrected because data is from simulation');
end

end

