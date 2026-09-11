function [midpoints,durations,integrationName] = tl_time_layers(model,count)
%TL_TIME_LAYERS One period, including every square-wave time interface.
T = model.modulation.period;
if strcmp(model.modulation.type,'square')
    phase = mod(model.modulation.phase,2*pi);
    boundaries = [0,2*pi*model.modulation.dutyCycle];
    times = mod(boundaries-phase,2*pi)/model.modulation.OmegaRadPerSec;
    times = times(times > 64*eps(T) & times < T-64*eps(T));
    times = unique(sort([0,times,T]));
    durations = diff(times);
    midpoints = times(1:end-1)+durations/2;
    integrationName = 'exact-square-layers';
else
    durations = repmat(T/count,1,count);
    midpoints = ((1:count)-0.5)*T/count;
    integrationName = 'midpoint-expm';
end
end
