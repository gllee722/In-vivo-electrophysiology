function [rateMat, binaryMat, tCenters] = binSpikePerTrial(spikeMat, binMs, startTimeMs, sampleStepMs, sigmaBins)
% spikeMat: trials x time, 0/1 (binary; if one spike spans multiple samples, it will be reduced to a single rising edge)
% binMs: bin size in ms (e.g., 10 ms)
% startTimeMs: start time of the window in ms (e.g., -1000)
% sampleStepMs: sampling interval in ms (must match the data; e.g., 1 ms)
% sigmaBins: standard deviation of Gaussian smoothing kernel in units of bins (e.g., 1)

if nargin < 5, sigmaBins = 1; end

[numTrials, numTimePoints] = size(spikeMat);

% --- 1) Remove duplicate rising edges (to ensure a single spike is counted once)
spikeMat = logical(spikeMat);
spikeMat(:,2:end) = spikeMat(:,2:end) & ~spikeMat(:,1:end-1);

% --- 2) Use only full bins (e.g., exactly 400 bins)
binSamples = round(binMs / sampleStepMs);       % number of samples per bin (e.g., 10)
nBins      = floor(numTimePoints / binSamples); % total number of full bins (e.g., 400)
usePts     = nBins * binSamples;                % number of samples used (e.g., 4000)
sp = spikeMat(:,1:usePts);

% --- 3) Compute spike counts per bin (vectorized, no loops)
sp = reshape(sp, numTrials, binSamples, nBins); % [trials x samples-per-bin x bins]
counts = squeeze(sum(sp, 2));                   % [trials x bins]

% --- 4) Convert counts to firing rate (Hz)
binSec  = (binSamples * sampleStepMs) / 1000;   % bin duration in seconds (e.g., 10 ms → 0.01 s)
rateMat = counts / binSec;                      % [trials x bins], firing rate in Hz

% --- 5) Binary representation: whether each bin contains any spike
binaryMat = counts > 0;                         % [trials x bins], 0/1

% --- 6) Gaussian smoothing (kernel normalized so mean rate is preserved)
if sigmaBins > 0
    halfW = max(1, round(5*sigmaBins));         % truncate kernel at ±5σ
    x = -halfW:halfW;
    g = exp(-(x.^2)/(2*sigmaBins^2));
    g = g / sum(g);                              % normalize kernel
    % Convolve along the time dimension for each trial
    rateMat = conv2(rateMat, g, 'same');         % preserves size [trials x bins]
end

% --- 7) Bin-center times in ms
t0 = startTimeMs;
tCenters = t0 + ( (0:nBins-1) + 0.5 ) * binSamples * sampleStepMs;
end