function [lam, err_gm] = PT_optimize(ytr, pmax, pmax2, Fs)
% -------------------------------------------------------------------------
% PT_OPTIMIZE
% Optimize pole-tracking parameters using GA.
%
% This function runs the GA optimizer using the PT cost function
% (through GAPT) and returns the optimal parameter vector.
%
% INPUTS
%   ytr   : training signals [M x N x trials]
%   pmax  : number of oscillatory pole pairs
%   pmax2 : number of real poles
%   Fs    : sampling frequency
%
% OUTPUTS
%   lam    : optimized parameter vector
%   err_gm : final GA cost value
%
% Parameter vector structure:
%
% lam =
% [ R2
%   R1
%   oscillatory radii (pmax)
%   real poles (pmax2)
%   oscillatory angles (pmax)
%   P0 ]
%
% Total parameters = 3 + 2*pmax + pmax2
% -------------------------------------------------------------------------

metric = 1;
ignore = 5;

%% ================= GA SETTINGS =================

fminconOptions = optimoptions(@fmincon,...
    'Display','iter',...
    'UseParallel',true,...
    'Algorithm','active-set');

ga_opts = gaoptimset( ...
    'TolFun',1e-12,...
    'StallGenLimit',20,...
    'Generations',80,...
    'Display','iter',...
    'UseParallel',true,...
    'HybridFcn',{@fmincon,fminconOptions});

%% ================= PARAMETER VECTOR =================

nvars = 3 + 2*pmax + pmax2;

% Lower bounds
LB = [ ...
    1e-5 ...                     % R2
    1e-5 ...                     % R1
    1e-5*ones(1,pmax) ...       % oscillatory radii
    -1*ones(1,pmax2) ...        % real poles
    zeros(1,pmax) ...           % oscillatory angles
    1e-5];                      % P0

% Upper bounds
UB = [ ...
    inf ...
    1 ...
    ones(1,pmax) ...
    ones(1,pmax2) ...
    pi*ones(1,pmax) ...
    inf];

%% ================= COST FUNCTION =================

h = @(X) GAPT(X, ytr, metric, ignore, pmax, pmax2, Fs);

%% ================= RUN GA =================

fprintf('Running PT optimization (GA)...\n');

[lam, err_gm] = ga(h, nvars, [], [], [], [], LB, UB, [], [], ga_opts);

fprintf('Optimization finished.\n');

end