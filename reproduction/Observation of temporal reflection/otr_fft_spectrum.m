function [f,X,info] = otr_fft_spectrum(t,x,varargin)
%OTR_FFT_SPECTRUM Toolbox-free, correctly scaled FFT spectrum.
%   [F,X,INFO] = OTR_FFT_SPECTRUM(T,X) transforms along the first
%   non-singleton dimension of X. T is uniformly spaced in seconds and F
%   is in Hz. Real input defaults to one-sided output; complex input to a
%   centred two-sided spectrum so negative-frequency reflection is kept.
%
%   Options (name/value):
%     'Dimension'  transform dimension (default first non-singleton)
%     'NFFT'       FFT length (default signal length)
%     'Window'     'none' (default), 'hann', or a numeric vector
%     'Detrend'    true/false (default false; removes mean)
%     'Spectrum'   'auto', 'onesided', 'twosided', or 'centered'
%
%   X is the complex amplitude spectrum (not a power spectrum). Scaling by
%   the window coherent gain makes a bin-centred sinusoid retain its peak
%   amplitude. INFO.magnitude, phase, power and frequency resolution are
%   included for plotting/retrieval.  无需 Signal Processing Toolbox。

if ~(isnumeric(t) && isreal(t) && isvector(t) && numel(t)>=2 && all(isfinite(t(:))))
    error('otr:fft:Time','t must be a finite real vector with at least two samples.');
end
if ~(isnumeric(x) && ~isempty(x) && all(isfinite(x(:))))
    error('otr:fft:Signal','x must be a finite nonempty numeric array.');
end
nT = numel(t);
dim = find(size(x)>1,1,'first');
if isempty(dim), dim = 1; end
nfft = nT; window = 'none'; detrendFlag = false; spectrum = 'auto';
if mod(numel(varargin),2) ~= 0
    error('otr:fft:Options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(char(varargin{k})); value = varargin{k+1};
    switch name
        case 'dimension', dim = value;
        case 'nfft', nfft = value;
        case 'window', window = value;
        case 'detrend', detrendFlag = value;
        case 'spectrum', spectrum = lower(char(value));
        otherwise, error('otr:fft:Options','Unknown option ''%s''.',name);
    end
end
if ~(isnumeric(dim) && isscalar(dim) && dim==round(dim) && dim>=1 && dim<=ndims(x))
    error('otr:fft:Dimension','Dimension must be a valid positive integer.');
end
if size(x,dim) ~= nT
    error('otr:fft:Length','numel(t) must equal size(x,Dimension).');
end
if ~(isnumeric(nfft) && isscalar(nfft) && isfinite(nfft) && nfft==round(nfft) && nfft>=nT)
    error('otr:fft:NFFT','NFFT must be an integer at least numel(t).');
end
if ~(islogical(detrendFlag) || (isnumeric(detrendFlag) && isscalar(detrendFlag)))
    error('otr:fft:Detrend','Detrend must be a logical scalar.');
end

tv = t(:);
dtv = diff(tv);
dt = mean(dtv);
if dt <= 0 || max(abs(dtv-dt)) > max(1e-12*abs(dt),64*eps(max(abs(tv))))
    error('otr:fft:Sampling','t must be strictly increasing and uniformly spaced.');
end

w = localWindow(window,nT);
gain = sum(w);
if abs(gain) < eps
    error('otr:fft:Window','Window has zero coherent gain.');
end
shape = ones(1,ndims(x)); shape(dim) = nT;
w = reshape(w,shape);
if logical(detrendFlag), x = x-mean(x,dim); end
y = fft(x.*w,nfft,dim)/gain;

if strcmp(spectrum,'auto')
    if isreal(x), spectrum = 'onesided'; else, spectrum = 'centered'; end
end
fs = 1/dt;
switch spectrum
    case 'onesided'
        if ~isreal(x)
            error('otr:fft:OneSidedComplex','One-sided scaling is only defined here for real signals.');
        end
        nKeep = floor(nfft/2)+1;
        ids = repmat({':'},1,ndims(y)); ids{dim} = 1:nKeep;
        X = y(ids{:});
        factor = ones(nKeep,1);
        if rem(nfft,2)==0, factor(2:end-1)=2; else, factor(2:end)=2; end
        shapeF = ones(1,ndims(X)); shapeF(dim)=nKeep;
        X = X.*reshape(factor,shapeF);
        f = (0:nKeep-1).'*fs/nfft;
    case 'twosided'
        X = y;
        f = (0:nfft-1).'*fs/nfft;
    case 'centered'
        X = fftshift(y,dim);
        f = ((0:nfft-1)-floor(nfft/2)).'*fs/nfft;
    otherwise
        error('otr:fft:Spectrum','Spectrum must be auto, onesided, twosided or centered.');
end

info = struct('magnitude',abs(X),'phase',angle(X),'power',abs(X).^2, ...
    'sampleInterval',dt,'sampleRate',fs,'frequencyResolution',fs/nfft, ...
    'nfft',nfft,'dimension',dim,'spectrum',spectrum,'coherentGain',gain/nT);
end

function w = localWindow(window,n)
if ischar(window) || (isstring(window) && isscalar(window))
    switch lower(char(window))
        case {'none','rectangular','rect'}
            w = ones(n,1);
        case {'hann','hanning'}
            if n==1, w=1; else, w=0.5-0.5*cos(2*pi*(0:n-1).'/(n-1)); end
        otherwise
            error('otr:fft:Window','Unknown window ''%s''.',char(window));
    end
elseif isnumeric(window) && isvector(window) && numel(window)==n && all(isfinite(window(:)))
    w = window(:);
else
    error('otr:fft:Window','Window must be none, hann, or a finite length-N vector.');
end
end
