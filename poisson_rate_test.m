function [p, RR, Cpre, Cpost, Tpre, Tpost] = poisson_rate_test(preCounts, postCounts, win_sec)
% Exact Poisson rate test (conditional binomial test)
%
% INPUTS:
%   preCounts  - vector, spike counts in the pre window for each trial
%   postCounts - vector, spike counts in the post window for each trial
%   win_sec    - scalar, duration of each window in seconds
%
% OUTPUTS:
%   p     - exact two-tailed p-value
%   RR    - rate ratio (post_rate / pre_rate)
%   Cpre  - total spike count in pre window
%   Cpost - total spike count in post window
%   Tpre  - total duration of pre window (s)
%   Tpost - total duration of post window (s)

% === Compute total spike counts ===
Cpre  = sum(preCounts);
Cpost = sum(postCounts);
N     = numel(preCounts);

% === Compute total time ===
Tpre  = N * win_sec;
Tpost = N * win_sec;

% === Exact conditional binomial test ===
Ctot = Cpre + Cpost;
p0   = Tpost / (Tpre + Tpost);

p_right = 1 - binocdf(Cpost-1, Ctot, p0);
p_left  = binocdf(Cpost, Ctot, p0);
p = 2 * min(p_left, p_right);   % two-tailed p-value

% === Rate ratio ===
RR = (Cpost/Tpost) / (Cpre/Tpre);
end