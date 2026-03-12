function [J, JJ] = GAPT(X, yy, metric, ignore, pmax, pmax2, Fs)
% -------------------------------------------------------------------------
% GAPT
% Applies the pole tracking cost function PT to all signals in a dataset.
%
% INPUT
%   X      : parameter vector for PT
%   yy     : signals [M x N x K]
%            M = number of channels/signals
%            N = samples
%            K = trials/segments
%   metric : cost metric flag (passed to PT)
%   ignore : number of initial samples ignored in cost
%   pmax   : number of oscillatory pole pairs
%   pmax2  : number of real poles
%   Fs     : sampling frequency
%
% OUTPUT
%   J  : mean cost across all signals and trials
%   JJ : matrix of costs per signal/trial [M x K]
% -------------------------------------------------------------------------

[M, ~, K] = size(yy);

% Preallocate
J = zeros(M, K);

% Loop over trials
for j = 1:K
    
    % Loop over signals
    for m = 1:M
        
        y = squeeze(yy(m,:,j));
        
        J(m,j) = PT(X, y, ignore, pmax, pmax2);
        
    end
    
end

% Global mean cost
J = mean(J(:));

end