function [poles, thmo_all] = extractpoles(lam, ignore, y, Fs, pmax, pmax2)
% -------------------------------------------------------------------------
% PT_EXTRACTPOLES
% Run SIM_PT on all signals using an already optimized parameter vector lam
% and extract pole features.
%
% INPUTS
%   lam    : optimized parameter vector from PT_optimize
%   y      : signals for pole extraction [M x T]
%   Fs     : sampling frequency
%   pmax   : number of oscillatory pole pairs
%   pmax2  : number of real poles
%
% OUTPUTS
%   poles   : concatenated pole feature matrix
%   thmo_all : cell array with raw tracked parameter trajectories per signal
%
% poles row layout:
%   [oscillatory radii    ]  -> pmax*M rows
%   [oscillatory freq     ]  -> pmax*M rows
%   [real poles           ]  -> pmax2*M rows
%   [covariance norm ptr  ]  -> M rows
% -------------------------------------------------------------------------


M = size(y,1);
T = size(y,2);

% Preallocate feature storage
polf  = cell(1,pmax);    % oscillatory radii
pola  = cell(1,pmax);    % oscillatory frequencies
polar = cell(1,pmax2);   % real poles

for k = 1:pmax
    polf{k} = nan(M,T);
    pola{k} = nan(M,T);
end

for k = 1:pmax2
    polar{k} = nan(M,T);
end

ptr = nan(M,T);

% Optional raw output
thmo_all = cell(M,1);

%% ================= RUN TRACKER =================
for m = 1:M
    yy = double(squeeze(y(m,:)));

    [thmo, ~, ~, ~, ptrr] = SIM_PT(lam, yy, ignore, pmax, pmax2, Fs);
    thmo_all{m} = thmo;

    % Oscillatory poles
    for k = 1:pmax
        polf{k}(m,:) = thmo(k,:);
        pola{k}(m,:) = abs(thmo(pmax + pmax2 + k,:)) * Fs / (2*pi);
    end

    % Real poles
    for k = 1:pmax2
        polar{k}(m,:) = thmo(pmax + k,:);
    end

    ptr(m,:) = ptrr(:).';
end

%% ================= CONCATENATE FEATURES =================
poles = [];

for k = 1:pmax
    poles = cat(1, poles, polf{k});
end

for k = 1:pmax
    poles = cat(1, poles, pola{k});
end

for k = 1:pmax2
    poles = cat(1, poles, polar{k});
end

poles = cat(1, poles, ptr);

end