function [y, x_clean, f_true_all, A_true_all, t, info] = ...
    generate_signals(M, Fs, Tsec, snr_db, params)

%% ---------------- CHECKS ----------------
nSin = params.nSin;

if ~isfield(params,'f_min_global')
    params.f_min_global = 5;
end

if ~isfield(params,'f_max_global')
    params.f_max_global = Fs/2 - 1;
end

%% ---------------- TIME ----------------
t = 0:1/Fs:Tsec-1/Fs;
N = length(t);

%% ---------------- PREALLOCATE ----------------
f_true_all = zeros(M,N);
A_true_all = zeros(M,N);
x_clean    = zeros(M,N);
y          = zeros(M,N);

info.mode = cell(M,1);
info.fmin = zeros(M,1);
info.fmax = zeros(M,1);

%% ---------------- GENERATE SIGNALS ----------------
for m = 1:M

    % Reuse parameter sets cyclically if M > nSin
    idx = mod(m-1, nSin) + 1;

    fmin = params.f_min(idx);
    fmax = params.f_max(idx);

    Amin = params.A_min(idx);
    Amax = params.A_max(idx);

    % Enforce global bounds
    fmin = max(fmin, params.f_min_global);
    fmax = min(fmax, params.f_max_global);
    if fmax <= fmin
        fmax = fmin + eps;
    end

    %% Choose mode
    if strcmpi(params.mode,'mixed')
        modes = {'linear','sinusoidal','abrupt'};
        this_mode = modes{randi(3)};
    else
        this_mode = params.mode;
    end

    info.mode{m} = this_mode;
    info.fmin(m) = fmin;
    info.fmax(m) = fmax;

    %% Generate trajectories
    [f_true, A_true] = local_generate_trajectory( ...
        t, Fs, fmin, fmax, Amin, Amax, this_mode, params);

    %% Enforce global frequency limits
    f_true = max(f_true, params.f_min_global);
    f_true = min(f_true, params.f_max_global);

    %% Build oscillator
    phi = cumsum(2*pi*f_true/Fs);
    x_m = A_true .* sin(phi);

    %% Store clean signal and truth
    x_clean(m,:)    = x_m;
    f_true_all(m,:) = f_true;
    A_true_all(m,:) = A_true;

    %% Add noise
    sigPow   = var(x_m);
    noisePow = sigPow/(10^(snr_db/10));
    y(m,:)   = x_m + sqrt(noisePow)*randn(1,N);

end

end


function [f_true, A_true] = local_generate_trajectory( ...
    t, Fs, fmin, fmax, Amin, Amax, this_mode, params)

N = length(t);

switch lower(this_mode)

    case 'linear'
        f_start = fmin + (fmax-fmin)*rand;
        f_end   = fmin + (fmax-fmin)*rand;
        f_true  = linspace(f_start,f_end,N);

        A_start = Amin + (Amax-Amin)*rand;
        A_end   = Amin + (Amax-Amin)*rand;
        A_true  = linspace(A_start,A_end,N);

    case 'sinusoidal'
        f_center = (fmin+fmax)/2;
        f_amp    = (fmax-fmin)/4;

        f_modHz = params.f_mod_min + ...
                  (params.f_mod_max-params.f_mod_min)*rand;

        f_true = f_center + ...
                 f_amp*sin(2*pi*f_modHz*t + 2*pi*rand);

        A_center = (Amin+Amax)/2;
        A_amp    = (Amax-Amin)/4;

        A_modHz = params.A_mod_min + ...
                  (params.A_mod_max-params.A_mod_min)*rand;

        A_true = A_center + ...
                 A_amp*sin(2*pi*A_modHz*t + 2*pi*rand);

    case 'abrupt'
        nSeg = randi([params.nSeg_min params.nSeg_max]);
        bp = round(linspace(1,N+1,nSeg+1));

        f_true = zeros(1,N);
        A_true = zeros(1,N);

        for s = 1:nSeg
            idx = bp(s):bp(s+1)-1;

            f_val = fmin + (fmax-fmin)*rand;
            A_val = Amin + (Amax-Amin)*rand;

            f_true(idx) = f_val;
            A_true(idx) = A_val;
        end

        win = max(3,round(params.abrupt_smooth_sec*Fs));
        f_true = movmean(f_true,win);
        A_true = movmean(A_true,win);

    otherwise
        error('Unknown mode: %s', this_mode);
end

end