function [spectrogramValue, f, tCenter] = stm_stft(signal, dt, ...
    windowLength, hop, nFFT)
%STM_STFT Toolbox-free short-time Fourier transform.

signal = signal(:);
if nargin < 3 || isempty(windowLength)
    windowLength = 128;
end
if nargin < 4 || isempty(hop)
    hop = floor(windowLength/8);
end
if nargin < 5 || isempty(nFFT)
    nFFT = 2^nextpow2(2*windowLength);
end

window = 0.5 - 0.5*cos(2*pi*(0:windowLength-1).' ...
    /(windowLength-1));
starts = 1:hop:(numel(signal)-windowLength+1);
nPositive = floor(nFFT/2)+1;
spectrogramValue = complex(zeros(nPositive,numel(starts)));

for q = 1:numel(starts)
    ids = starts(q):(starts(q)+windowLength-1);
    localSpectrum = fft(signal(ids).*window,nFFT);
    spectrogramValue(:,q) = localSpectrum(1:nPositive);
end

f = (0:nPositive-1).'/(nFFT*dt);
tCenter = ((starts-1)+(windowLength-1)/2)*dt;
end
