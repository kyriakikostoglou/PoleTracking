function [thmo, thm, e, J, ptr] = SIM_PT(X, Y, ignore, pmax, pmax2, Fs)
% -------------------------------------------------------------------------
% SIM_PT
% Simulation / tracking version of PT.
%
% INPUTS
%   X      : parameter vector
%   Y      : signal row vector
%   ignore : number of initial samples to ignore in cost
%   pmax   : number of oscillatory pole pairs
%   pmax2  : number of real poles
%   Fs     : sampling frequency
%
% OUTPUTS
%   thmo : tracked parameters (optionally ordered by frequency)
%   thm  : parameter trajectories
%   e    : prediction error
%   J    : normalized prediction error
%   ptr  : covariance norm over time
% -------------------------------------------------------------------------

pim = pmax;          % oscillatory pairs
pr  = pmax2;         % real poles
p   = pim + pr;

% Parameter vector:
% X = [R2, R1, osc_radii(1:pim), real_poles(1:pr), osc_angles(1:pim), P0]
th = X(3:2 + 2*pim + pr).';
totPar = length(th);

R2 = X(1);
R1 = X(2) * eye(totPar);
P0 = X(end);
P  = P0 * eye(totPar);

N = size(Y,2);

% Preallocate
e    = zeros(1,N);
yy1  = zeros(pim,N);
yy2  = zeros(pim,N);
yyr  = zeros(pr,N);
epi  = zeros(p,N);
phi  = zeros(totPar,1);
ptr  = zeros(N,1);

thm  = zeros(totPar,N);
thmo = complex(nan(totPar,N));   % same size as thm

for k = 3:N

    % -------------------------------------------------------------
    % Auxiliary states for oscillatory poles
    % -------------------------------------------------------------
    for i = 1:pim
        a1 = -2 * th(i) * cos(th(i+p));
        a2 = th(i)^2;

        yy1(i,k) = -a1 * yy1(i,k-1) - a2 * yy1(i,k-2) + e(k-1);
        yy2(i,k) = -a1 * yy2(i,k-1) - a2 * yy2(i,k-2) + e(k-2);
    end

    % -------------------------------------------------------------
    % Auxiliary states for real poles
    % -------------------------------------------------------------
    for i = 1:pr
        a1 = -th(i+pim);
        yyr(i,k) = -a1 * yyr(i,k-1) + e(k-1);
    end

    % -------------------------------------------------------------
    % Build regressor phi
    % -------------------------------------------------------------
    oscCount = 0;
    for i = 1:totPar
        if i <= pim
            % derivative wrt oscillatory radii
            phi(i) = -(-2*cos(th(i+p))*yy1(i,k) + 2*th(i)*yy2(i,k));
        elseif i <= p
            % derivative wrt real poles
            phi(i) = yyr(i-pim,k);
        else
            % derivative wrt oscillatory angles
            oscCount = oscCount + 1;
            phi(i) = -(2*sin(th(i))*yy1(oscCount,k));
        end
    end

    % -------------------------------------------------------------
    % Cascaded prediction error
    % -------------------------------------------------------------
    if pim > 0
        count = 0;

        % Oscillatory sections
        for i = 1:pim
            count = count + 1;
            a1 = -2 * th(i) * cos(th(i+p));
            a2 = th(i)^2;

            if i == 1
                epi(count,k) = Y(k) + a1*Y(k-1) + a2*Y(k-2);
            else
                epi(count,k) = epi(count-1,k) ...
                             + a1*epi(count-1,k-1) ...
                             + a2*epi(count-1,k-2);
            end
        end

        % Real-pole sections
        for i = pim+1:pim+pr
            count = count + 1;
            a1 = -th(i);
            epi(count,k) = epi(count-1,k) + a1*epi(count-1,k-1);
        end

    else
        count = 0;

        for i = 1:pr
            count = count + 1;
            a1 = -th(i);

            if i == 1
                epi(count,k) = Y(k) + a1*Y(k-1);
            else
                epi(count,k) = epi(count-1,k) + a1*epi(count-1,k-1);
            end
        end
    end

    e(k) = epi(count,k);

    % -------------------------------------------------------------
    % Kalman update
    % -------------------------------------------------------------
    P = P + R1;

    pphi  = P * phi;
    denom = R2 + phi' * P * phi;
    K     = pphi / denom;

    P = (eye(size(P)) - K * phi') * P;
    ptr(k) = norm(P,'fro');

    th = th + K * e(k);

    % -------------------------------------------------------------
    % Stability / boundedness corrections
    % -------------------------------------------------------------
    if pim > 0
        % keep oscillatory radii inside unit circle
        if any(abs(th(1:pim)) > 1)
            inds = find(abs(th(1:pim)) > 1);
            temp = 1 ./ conj(th(inds) .* exp(1i * th(inds+p)));
            th(inds)   = abs(temp);
            th(inds+p) = abs(angle(temp));
        end

        % wrap angles to [-pi, pi]
        if any(abs(th(p+1:end)) > pi)
            inds = find(abs(th(p+1:end)) > pi);
            temp = th(p + inds);
            th(p + inds) = abs(wrapToPi(temp));
        end

        % keep real poles inside unit circle
        if pr > 0
            realPart = th(pim+1:p);
            if any(abs(realPart) > 1)
                inds = find(abs(realPart) > 1);
                th(pim + inds) = 1 ./ realPart(inds);
            end
        end

        th(1:pim) = abs(th(1:pim));

    else
        if any(abs(th(1:p)) > 1)
            inds = find(abs(th(1:p)) > 1);
            th(inds) = 1 ./ th(inds);
        end
    end

    % -------------------------------------------------------------
    % Optional ordering by oscillatory frequency
    % -------------------------------------------------------------
    if pim > 1
        [~, idx] = sort(abs(th(p+1:end)) * Fs / (2*pi), 'ascend');

        temp = th;
        temp(1:pim)   = th(idx);
        temp(p+1:end) = th(p + idx);

        thm(:,k) = temp;
    else
        thm(:,k) = th;
    end
end

thmo = thm;
J = norm(e(ignore:end)) / norm(Y(ignore:end));
end