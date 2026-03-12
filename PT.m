function [J, Up] = PT(X, Y, ignore, pmax, pmax2)
% -------------------------------------------------------------------------
% PT
% Kalman filter based pole tracking cost function
%
% INPUTS:
%   X      : parameter vector
%   Y      : signal row vector
%   ignore : number of initial samples to ignore in cost
%   pmax   : number of oscillatory pole pairs
%   pmax2  : number of real poles
%
% OUTPUTS:
%   J      : normalized prediction error
%   Up     : prediction error after ignored samples
% -------------------------------------------------------------------------

pim = pmax;      % number of oscillatory poles
pr  = pmax2;     % number of real poles
p   = pim + pr;

% Parameter vector:
% X = [R2, R1, oscillatory radii(1:pim), real poles(1:pr), angles(1:pim), P0]
th  = X(3:2 + 2*pim + pr).';
totPar = length(th);

R2  = X(1);
R1  = X(2) * eye(totPar);
P0  = X(end);
P   = P0 * eye(totPar);

N = size(Y, 2);

% Preallocation
e   = zeros(1, N);          % innovation / prediction error
yy1 = zeros(pim, N);        % derivative-related auxiliary states
yy2 = zeros(pim, N);
yyr = zeros(pr, N);         % auxiliary states for real poles
phi = zeros(totPar, 1);     % regressor
epi = zeros(p, N);          % cascaded prediction errors
thm = zeros(totPar, N);     % tracked parameters over time

for k = 3:N

    % -------------------------------------------------------------
    % Build auxiliary states for oscillatory and real poles
    % -------------------------------------------------------------
    for i = 1:pim
        a1 = -2 * th(i) * cos(th(i + p));
        a2 = th(i)^2;

        yy1(i,k) = -a1 * yy1(i,k-1) - a2 * yy1(i,k-2) + e(k-1);
        yy2(i,k) = -a1 * yy2(i,k-1) - a2 * yy2(i,k-2) + e(k-2);
    end

    for i = 1:pr
        a1 = -th(i + pim);
        yyr(i,k) = -a1 * yyr(i,k-1) + e(k-1);
    end

    % -------------------------------------------------------------
    % Build regressor phi
    % First pim entries: derivatives wrt oscillatory radii
    % Next pr entries  : derivatives wrt real poles
    % Last pim entries : derivatives wrt oscillatory angles
    % -------------------------------------------------------------
    oscCount = 0;
    for i = 1:totPar
        if i <= pim
            phi(i) = -(-2*cos(th(i+p))*yy1(i,k) + 2*th(i)*yy2(i,k));
        elseif i <= p
            phi(i) = yyr(i-pim,k);
        else
            oscCount = oscCount + 1;
            phi(i) = -(2*sin(th(i))*yy1(oscCount,k));
        end
    end

    % -------------------------------------------------------------
    % Build cascaded error epi
    % -------------------------------------------------------------
    if pim > 0
        count = 0;

        % oscillatory sections
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

        % real-pole sections
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

    pphi = P * phi;
    denom = R2 + phi' * P * phi;
    K = pphi / denom;

    P = (eye(size(P)) - K * phi') * P;
    th = th + K * e(k);

    % -------------------------------------------------------------
    % Stability / boundedness corrections
    % -------------------------------------------------------------
    if pim > 0
        % oscillatory radii should remain <= 1
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

        % real poles should remain inside unit circle
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

    thm(:,k) = th;
end

Up = e(ignore:end);
J  = norm(Up) / norm(Y(ignore:end));
end