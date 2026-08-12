function schedule = otr_switch_schedule(t, specification)
%OTR_SWITCH_SCHEDULE Create switch weights and exact event metadata.
%   SCHEDULE = OTR_SWITCH_SCHEDULE(T,SPECIFICATION) samples a spatially
%   uniform (or per-switch) switching law on the strictly increasing time
%   grid T.  Event times must lie on T to roundoff accuracy; this avoids
%   silently applying a temporal boundary between integration steps.
%
%   SPECIFICATION can be a scalar structure with fields
%       eventTimes    increasing vector of event times
%       states        scalar/uniform sequence with numel(eventTimes)+1
%                     entries, or an nSwitch-by-(nEvents+1) matrix
%       nSwitch       number of switches (default 1)
%       riseTime      nonnegative scalar; zero gives exact steps
%       transition    'linear' (default) or 'smoothstep'
%       snapTolerance optional absolute grid tolerance
%
%   Convenience fields initialState and eventStates may replace states.
%   A string mode can also be supplied in SPECIFICATION.mode:
%       'off-on'      one event, 0 -> 1
%       'on-off'      one event, 1 -> 0
%       'slab'        two events, 0 -> 1 -> 0
%       'inverse-slab'two events, 1 -> 0 -> 1
%
%   For finite riseTime, transitions start at eventTimes and finish after
%   riseTime.  Events may not overlap.  Exact event indices are returned in
%   SCHEDULE.eventIndices.

if nargin < 2 || isempty(specification)
    specification = struct();
end
if ~(isnumeric(t) && isreal(t) && isvector(t) && numel(t) >= 2 && ...
        all(isfinite(t(:))))
    error('otr_switch_schedule:InvalidTime', ...
        'T must be a finite real vector containing at least two samples.');
end
t = double(t(:).');
if any(diff(t) <= 0)
    error('otr_switch_schedule:InvalidTime', ...
        'T must be strictly increasing.');
end
if ~isstruct(specification) || ~isscalar(specification)
    error('otr_switch_schedule:InvalidSpecification', ...
        'SPECIFICATION must be a scalar structure.');
end

nSwitch = local_field(specification,'nSwitch',1);
if ~(isnumeric(nSwitch) && isscalar(nSwitch) && isfinite(nSwitch) && ...
        nSwitch >= 1 && nSwitch == round(nSwitch))
    error('otr_switch_schedule:InvalidSwitchCount', ...
        'nSwitch must be a positive integer.');
end
eventTimesRequested = double(local_field(specification,'eventTimes',[]));
eventTimesRequested = eventTimesRequested(:).';
if any(~isfinite(eventTimesRequested)) || any(diff(eventTimesRequested) <= 0)
    error('otr_switch_schedule:InvalidEvents', ...
        'eventTimes must be finite and strictly increasing.');
end
if any(eventTimesRequested < t(1)) || any(eventTimesRequested > t(end))
    error('otr_switch_schedule:EventOutsideGrid', ...
        'Every event time must lie inside the supplied time grid.');
end
nEvents = numel(eventTimesRequested);

if isfield(specification,'states') && ~isempty(specification.states)
    states = double(specification.states);
elseif isfield(specification,'mode') && ~isempty(specification.mode)
    mode = lower(char(specification.mode));
    switch mode
        case 'off-on'
            states = [0 1];
        case 'on-off'
            states = [1 0];
        case 'slab'
            states = [0 1 0];
        case 'inverse-slab'
            states = [1 0 1];
        otherwise
            error('otr_switch_schedule:InvalidMode', ...
                'Unknown switching mode ''%s''.', mode);
    end
else
    initialState = local_field(specification,'initialState',0);
    eventStates = local_field(specification,'eventStates',ones(1,nEvents));
    states = [initialState(:), eventStates];
end

if isvector(states)
    if numel(states) ~= nEvents+1
        error('otr_switch_schedule:InvalidStates', ...
            'A state sequence must have numel(eventTimes)+1 entries.');
    end
    states = repmat(states(:).',nSwitch,1);
elseif ~isequal(size(states),[nSwitch,nEvents+1])
    error('otr_switch_schedule:InvalidStates', ...
        'states must be a vector or an nSwitch-by-(nEvents+1) matrix.');
end
if any(~isfinite(states(:))) || any(states(:)<0) || any(states(:)>1)
    error('otr_switch_schedule:InvalidStates', ...
        'All states/connection weights must lie in [0,1].');
end

dt = diff(t);
defaultTolerance = max(64*eps(max(1,max(abs(t)))), ...
    1e-8*min(dt));
snapTolerance = local_field(specification,'snapTolerance',defaultTolerance);
if ~(isnumeric(snapTolerance) && isscalar(snapTolerance) && ...
        isfinite(snapTolerance) && snapTolerance >= 0)
    error('otr_switch_schedule:InvalidTolerance', ...
        'snapTolerance must be a nonnegative finite scalar.');
end
eventIndices = zeros(1,nEvents);
eventTimes = zeros(1,nEvents);
for event = 1:nEvents
    [distance,index] = min(abs(t-eventTimesRequested(event)));
    if distance > snapTolerance
        error('otr_switch_schedule:EventNotOnGrid', ...
            ['Event %.16g is not on the integration grid (nearest sample ' ...
             '%.16g, error %.3g).'], eventTimesRequested(event), ...
             t(index), distance);
    end
    eventIndices(event) = index;
    eventTimes(event) = t(index);
end
if any(diff(eventIndices) <= 0)
    error('otr_switch_schedule:EventsShareSample', ...
        'Each event must map to a distinct increasing grid sample.');
end

riseTime = local_field(specification,'riseTime',0);
if ~(isnumeric(riseTime) && isscalar(riseTime) && isreal(riseTime) && ...
        isfinite(riseTime) && riseTime >= 0)
    error('otr_switch_schedule:InvalidRiseTime', ...
        'riseTime must be a nonnegative finite scalar.');
end
transition = lower(char(local_field(specification,'transition','linear')));
if ~ismember(transition,{'linear','smoothstep'})
    error('otr_switch_schedule:InvalidTransition', ...
        'transition must be ''linear'' or ''smoothstep''.');
end
if riseTime > 0 && nEvents > 1 && any(diff(eventTimes) < riseTime-snapTolerance)
    error('otr_switch_schedule:OverlappingTransitions', ...
        'Finite-rise transitions may not overlap.');
end

weights = repmat(states(:,1),1,numel(t));
for event = 1:nEvents
    before = states(:,event);
    after = states(:,event+1);
    if riseTime == 0
        weights(:,eventIndices(event):end) = repmat(after,1, ...
            numel(t)-eventIndices(event)+1);
    else
        transitionMask = t >= eventTimes(event) & ...
            t < eventTimes(event)+riseTime;
        fraction = (t(transitionMask)-eventTimes(event))/riseTime;
        if strcmp(transition,'smoothstep')
            fraction = fraction.^2 .* (3-2*fraction);
        end
        if any(transitionMask)
            weights(:,transitionMask) = before + ...
                (after-before)*fraction;
        end
        weights(:,t >= eventTimes(event)+riseTime) = repmat(after,1, ...
            nnz(t >= eventTimes(event)+riseTime));
    end
end

schedule = struct();
schedule.time = t;
schedule.weights = weights;
schedule.states = states;
schedule.eventTimesRequested = eventTimesRequested;
schedule.eventTimes = eventTimes;
schedule.eventIndices = eventIndices;
schedule.riseTime = riseTime;
schedule.transition = transition;
schedule.nSwitch = nSwitch;
schedule.isInstantaneous = riseTime == 0;
schedule.eventType = zeros(nSwitch,nEvents);
for event = 1:nEvents
    delta = states(:,event+1)-states(:,event);
    schedule.eventType(delta>0,event) = 1;
    schedule.eventType(delta<0,event) = -1;
end

end

function value = local_field(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
