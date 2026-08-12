function gated = otr_gate_signal(time, signal, limits, taperFraction)
%OTR_GATE_SIGNAL Extract, baseline-correct and softly gate a waveform.
%   G = OTR_GATE_SIGNAL(T, X, [T0 T1]) returns a struct with samples in the
%   requested closed interval.  A line connecting robust endpoint means is
%   removed, then a cosine taper is applied only near the two boundaries.
%   This transparent operation avoids hidden filtering and requires no
%   Signal Processing Toolbox.
%
%   G = OTR_GATE_SIGNAL(..., TAPERFRACTION) sets the tapered fraction at
%   each end (default 0.10, valid range 0 to 0.5).  Set it to zero to retain
%   a rectangular gate.  Fields include time, raw, baseline, detrended,
%   window, signal, limits and dt.

    if nargin < 3 || numel(limits) ~= 2
        error('otr:BadGateLimits', 'limits must contain [start end].');
    end
    if nargin < 4 || isempty(taperFraction)
        taperFraction = 0.10;
    end
    if taperFraction < 0 || taperFraction > 0.5
        error('otr:BadTaper', 'taperFraction must lie between 0 and 0.5.');
    end
    if ~isvector(time) || ~isvector(signal) || numel(time) ~= numel(signal)
        error('otr:BadWaveform', 'time and signal must be equal-length vectors.');
    end

    time = time(:);
    signal = signal(:);
    limits = sort(limits(:).');
    use = isfinite(time) & isfinite(signal) & time >= limits(1) & time <= limits(2);
    t = time(use);
    x = signal(use);
    if numel(t) < 4
        error('otr:EmptyGate', 'Gate [%.6g %.6g] contains fewer than 4 samples.', limits);
    end

    n = numel(x);
    endpointCount = max(3, min(floor(n / 4), round(0.08 * n)));
    startMean = local_trimmed_mean(x(1:endpointCount));
    endMean = local_trimmed_mean(x(end - endpointCount + 1:end));
    baseline = startMean + (endMean - startMean) * (0:n - 1).' / max(1, n - 1);
    detrended = x - baseline;

    window = ones(n, 1);
    edgeCount = min(floor(n / 2), round(taperFraction * n));
    if edgeCount > 0
        edge = 0.5 * (1 - cos(pi * (0:edgeCount - 1).' / edgeCount));
        window(1:edgeCount) = edge;
        window(end - edgeCount + 1:end) = flipud(edge);
    end

    dtValues = diff(t);
    dtValues = dtValues(isfinite(dtValues) & dtValues > 0);
    if isempty(dtValues)
        dt = NaN;
    else
        dt = median(dtValues);
    end

    gated = struct();
    gated.time = t;
    gated.raw = x;
    gated.baseline = baseline;
    gated.detrended = detrended;
    gated.window = window;
    gated.signal = detrended .* window;
    gated.limits = limits;
    gated.dt = dt;
end

function value = local_trimmed_mean(x)
    x = sort(x(isfinite(x)));
    if isempty(x)
        value = 0;
        return;
    end
    trim = floor(0.10 * numel(x));
    if 2 * trim < numel(x)
        x = x(trim + 1:end - trim);
    end
    value = mean(x);
end
