function [Rest, Ipred, prL1S, rest] = allFilSmoothGrow(Rgrid, m, eta, nday, p0, Lam, Iloc, distvals)

% Assumptions and notes
% - runs filter and smoother in one go
% - includes growth rates with confidence

% EpiFilter estimates using local cases
[~, RlowF, RhighF, RmeanF, pR, pRup, pstate] = runEpiFilterSm(Rgrid, m, eta, nday, p0, Lam, Iloc);

% EpiFilter one-step-ahead predictions 
[predF, predIntF] = recursPredictMax(Rgrid, pR, Lam, RmeanF, 2*max(Iloc));
predlowF = predIntF(:,1)'; predhighF = predIntF(:,2)';

% EpiSmoother estimates for single trajectory
[~, RlowS, RhighS, RmeanS, qR] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);

% EpiSmoother one-step-ahead predictions 
[predS, predIntS] = recursPredictMax(Rgrid, qR, Lam, RmeanS, 2*max(Iloc));
predlowS = predIntS(:,1)'; predhighS = predIntS(:,2)';

% For probabilities above or below 1
id1 = find(Rgrid <= 1, 1, 'last'); prL1S = zeros(1, nday); 
% Update prior to posterior sequentially
for i = 2:nday
    % Posterior CDF and prob R <= 1
    Rcdf = cumsum(qR(i, :)); prL1S(i) = Rcdf(id1);
end

% Growth rates from distributions over Rgrid if gamma SD
if distvals.type == 2
    % Gamma distribution shape and scale
    shape = distvals.pm;
    scale = distvals.omega/shape;
    
    % Growth rate from Wallinga-Lipsitch
    rmeanF = (RmeanF.^(1/shape) - 1)/scale;
    rmeanS = (RmeanS.^(1/shape) - 1)/scale;
    % Confidence intervals also direct as keep same grid point
    rlowF = (RlowF.^(1/shape) - 1)/scale;
    rlowS = (RlowS.^(1/shape) - 1)/scale;
    rhighF = (RhighF.^(1/shape) - 1)/scale;
    rhighS = (RhighS.^(1/shape) - 1)/scale;
    
    % Output data structures for growth estimates
    rest.mean = [rmeanF' rmeanS'];
    rest.low = [rlowF' rlowS'];
    rest.high = [rhighF' rhighS'];
else
    rest = []; disp('Gamma serial interval not used');
end

% Output data structures for R estimates
Rest.mean = [RmeanF' RmeanS'];
Rest.low = [RlowF' RlowS'];
Rest.high = [RhighF' RhighS'];

% Output data structures for I predictions
Ipred.mean = [predF' predS'];
Ipred.low = [predlowF' predlowS'];
Ipred.high = [predhighF' predhighS'];