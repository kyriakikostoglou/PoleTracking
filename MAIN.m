clc; close all; clear;

%% ================= SETTINGS =================
Fs     = 200;     % original generation sampling rate
Tsec   = 20;
snr_db = 20;
M      = 3;       % number of signals to generate

%% ================= SIGNAL PARAMETERS =================
% One sinusoid per signal, parameter sets reused if M > params.nSin
params.nSin = 3;

params.f_min = [8 14 22];
params.f_max = [12 18 28];

params.A_min = [0.4 0.2 0.15];
params.A_max = [1.0 0.6 0.4];

params.mode = 'mixed';   % 'linear' | 'sinusoidal' | 'abrupt' | 'mixed'

params.f_mod_min = 0.05;
params.f_mod_max = 0.20;
params.A_mod_min = 0.05;
params.A_mod_max = 0.20;

params.nSeg_min = 3;
params.nSeg_max = 6;
params.abrupt_smooth_sec = 0.05;

params.f_min_global = 5;
params.f_max_global = 40;

%% ================= GENERATE SIGNALS =================
[y, x_clean, f_true_all, A_true_all, t, info] = ...
    generate_signals(M, Fs, Tsec, snr_db, params);

%% ================= SHOW MODES =================
disp('Modes used per signal:')
for m = 1:M
    fprintf('Signal %d -> %s | band [%g %g] Hz\n', ...
        m, info.mode{m}, info.fmin(m), info.fmax(m));
end

%% ================= RESAMPLE TO 64 Hz (NO ANTI-ALIASING) =================
Fs_old = Fs;
Fs = 64;

y          = resample(y', Fs, Fs_old)';
x_clean    = resample(x_clean', Fs, Fs_old)';
f_true_all = resample(f_true_all', Fs, Fs_old)';
A_true_all = resample(A_true_all', Fs, Fs_old)';

N = size(y,2);
t = (0:N-1)/Fs;

%% ================= PT SETTINGS =================
pmax  = 1;
pmax2 = 0;

ignore = round(0.5 * Fs);
Ttrain = 10;
Ntrain = round(Ttrain * Fs);

%% ================= PREPARE TRAINING DATA =================
ytrain = y(:,1:Ntrain);
ytr    = reshape(ytrain,[M Ntrain 1]);

%% ================= OPTIMIZE PT =================
[lam, err_gm] = PT_optimize(ytr, pmax, pmax2, Fs);

fprintf('Final PT cost: %.6f\n', err_gm);
disp('Optimized parameters:')
disp(lam)

%% ================= RUN POLE TRACKING =================
[poles, thmo_all] = extractpoles(lam, ignore, y, Fs, pmax, pmax2);

%% ================= DECODE PT OUTPUT =================
offset = 0;

r_est = cell(1,pmax);
for k = 1:pmax
    r_est{k} = poles(offset + (1:M), :);
    offset = offset + M;
end

f_est = cell(1,pmax);
for k = 1:pmax
    f_est{k} = poles(offset + (1:M), :);
    offset = offset + M;
end

ptr = poles(offset + (1:M), :);

%% ================= PLOT SIGNALS =================
figure
for m = 1:M
    subplot(M,1,m)
    plot(t, y(m,:), 'LineWidth', 1); hold on
    plot(t, x_clean(m,:), 'LineWidth', 1.4)
    ylabel(sprintf('Sig %d', m))
    title(sprintf('Signal %d | mode: %s', m, info.mode{m}))
    xlim([t(1) t(end)])
end
xlabel('Time (s)')
sgtitle('Signals After Resampling to 64 Hz')
legend('Noisy','Clean')

%% ================= PLOT FREQUENCY TRACKING =================
t_plot = t(ignore:end);

figure
for m = 1:M
    subplot(M,1,m)
    hold on

    plot(t_plot, f_true_all(m,ignore:end), 'LineWidth', 2)

    for k = 1:pmax
        plot(t_plot, f_est{k}(m,ignore:end), '--', 'LineWidth', 1.4)
    end

    ylabel('Hz')
    title(sprintf('Signal %d', m))
    xlim([t_plot(1) t_plot(end)])
end
xlabel('Time (s)')
sgtitle('True vs PT Estimated Frequencies')
legend('True','Estimated')

%% ================= PLOT POLE MAGNITUDE =================
figure
for m = 1:M
    subplot(M,1,m)
    plot(t_plot, r_est{1}(m,ignore:end), 'LineWidth', 1.4)
    ylabel('Pole magnitude')
    title(sprintf('Signal %d', m))
    xlim([t_plot(1) t_plot(end)])
end
xlabel('Time (s)')
sgtitle('Pole Magnitude')

%% ================= PLOT COVARIANCE =================
figure
for m = 1:M
    subplot(M,1,m)
    plot(t_plot, ptr(m,ignore:end), 'LineWidth', 1.4)
    ylabel('||P||')
    title(sprintf('Signal %d', m))
    xlim([t_plot(1) t_plot(end)])
end
xlabel('Time (s)')
sgtitle('Covariance Norm')