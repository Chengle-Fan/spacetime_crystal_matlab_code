function result = otr_simulate_circuit(model, t, source, schedule, options)
%OTR_SIMULATE_CIRCUIT Integrate the discrete switched-line MNA equations.
%   RESULT = OTR_SIMULATE_CIRCUIT(MODEL,T,SOURCE,SCHEDULE,OPTIONS) solves
%
%       C(w) xdot + G(w) x = B*SOURCE(t)
%
%   with the trapezoidal rule.  MODEL is returned by OTR_BUILD_MNA.  SOURCE
%   is a function handle, a scalar, or one sample per entry of T.  SCHEDULE
%   is returned by OTR_SWITCH_SCHEDULE; a schedule specification structure
%   is also accepted and converted automatically.
%
%   The integration convention is deliberately explicit: an interval is
%   integrated with the state that exists immediately after its left end.
%   At its right end, a switch event is applied exactly on the grid.  Thus
%   RESULT.state(:,n) is x(t_n^+), while event records also retain x(t_n^-).
%
%   In ideal-capacitance mode, decreasing connection weights first remove
%   capacitors with voltage continuous.  Increasing weights then impose
%       C_before * v^- = C_after * v^+,
%   while every inductor current is continuous.  A finite rise is treated
%   as a sequence of these physically explicit incremental charge-sharing
%   events.  In RLC mode, weights couple a fixed internal RLC topology; an
%   OFF event isolates and clears the removed branch as in the experiment.
%
%   OPTIONS fields:
%       initialState       nState-vector (default all zero)
%       method             only 'trapezoidal' is currently supported
%       boundaryTolerance  relative validation tolerance (default 5e-10)
%       verifyBoundary     error when an exact boundary test fails (true)
%       storeState         retain the full state history (true)

if nargin < 5 || isempty(options)
    options = struct();
end
if ~isstruct(model) || ~isscalar(model) || ~isfield(model,'isMNA') || ...
        ~model.isMNA
    error('otr_simulate_circuit:InvalidModel', ...
        'MODEL must be returned by otr_build_mna.');
end
if ~(isnumeric(t) && isreal(t) && isvector(t) && numel(t)>=2 && ...
        all(isfinite(t(:))))
    error('otr_simulate_circuit:InvalidTime', ...
        'T must be a finite real vector with at least two entries.');
end
t = double(t(:).');
if any(diff(t)<=0)
    error('otr_simulate_circuit:InvalidTime', ...
        'T must be strictly increasing.');
end
if ~isstruct(options) || ~isscalar(options)
    error('otr_simulate_circuit:InvalidOptions', ...
        'OPTIONS must be a scalar structure.');
end
method = lower(char(local_field(options,'method','trapezoidal')));
if ~strcmp(method,'trapezoidal')
    error('otr_simulate_circuit:InvalidMethod', ...
        'Only the trapezoidal method is supported.');
end
boundaryTolerance = local_field(options,'boundaryTolerance',5e-10);
if ~(isnumeric(boundaryTolerance) && isscalar(boundaryTolerance) && ...
        isfinite(boundaryTolerance) && boundaryTolerance>0)
    error('otr_simulate_circuit:InvalidTolerance', ...
        'boundaryTolerance must be a positive finite scalar.');
end
verifyBoundary = local_logical(options,'verifyBoundary',true);
storeState = local_logical(options,'storeState',true);

nTime = numel(t);
sourceValues = local_source_values(source,t);

if nargin < 4 || isempty(schedule)
    schedule = otr_switch_schedule(t,struct('nSwitch',model.nSwitches));
elseif ~isfield(schedule,'weights') || ~isfield(schedule,'time')
    schedule.nSwitch = model.nSwitches;
    schedule = otr_switch_schedule(t,schedule);
end
schedule = local_validate_schedule(schedule,t,model.nSwitches);
weights = schedule.weights;

x0 = local_field(options,'initialState',zeros(model.nState,1));
if ~(isnumeric(x0) && isreal(x0) && isvector(x0) && ...
        numel(x0)==model.nState && all(isfinite(x0(:))))
    error('otr_simulate_circuit:InvalidInitialState', ...
        'initialState must be a finite real vector of length %d.',model.nState);
end
x = zeros(model.nState,nTime);
x(:,1) = double(x0(:));

eventMask = false(1,nTime);
if isfield(schedule,'eventIndices')
    eventMask(schedule.eventIndices) = true;
end
% A finite ramp has a small boundary update at every changing sample.
changeMask = [false, any(abs(diff(weights,1,2)) > ...
    64*eps(max(1,max(abs(weights(:))))),1)];
boundaryMask = eventMask | changeMask;
eventRecords = repmat(local_empty_event(model.nState,model.nSwitches), ...
    1,nnz(boundaryMask));
recordIndex = 0;

for sample = 2:nTime
    h = t(sample)-t(sample-1);
    oldWeight = weights(:,sample-1);
    newWeight = weights(:,sample);
    [Cleft,Gleft] = local_matrices(model,oldWeight);
    left = Cleft/h + Gleft/2;
    rhs = (Cleft/h-Gleft/2)*x(:,sample-1) + ...
        model.B*((sourceValues(sample-1)+sourceValues(sample))/2);
    xMinus = left \ rhs;
    if any(~isfinite(xMinus))
        error('otr_simulate_circuit:IntegrationFailure', ...
            'The MNA solve produced non-finite state values at sample %d.',sample);
    end

    if boundaryMask(sample)
        [xPlus,boundary] = local_apply_boundary(model,xMinus, ...
            oldWeight,newWeight);
        recordIndex = recordIndex+1;
        eventRecords(recordIndex) = local_make_event(t(sample),sample, ...
            oldWeight,newWeight,xMinus,xPlus,boundary);
        if verifyBoundary && boundary.relativeResidual > boundaryTolerance
            error('otr_simulate_circuit:BoundaryViolation', ...
                ['Temporal boundary residual %.3g exceeds tolerance %.3g ' ...
                 'at t=%.16g.'],boundary.relativeResidual, ...
                 boundaryTolerance,t(sample));
        end
        x(:,sample) = xPlus;
    else
        x(:,sample) = xMinus;
    end
end

v = x(model.index.voltage,:);
iLine = x(model.index.lineCurrent,:);
result = struct();
result.time = t;
result.source = sourceValues;
result.weights = weights;
result.voltage = v;
result.lineCurrent = iLine;
result.inputVoltage = v(1,:);
result.outputVoltage = v(end,:);
result.eventRecords = eventRecords;
result.schedule = schedule;
result.method = method;
result.modelPhysicsMode = model.physicsMode;
result.boundaryConvention = 'state(:,n)=x(t_n^+)';
if strcmp(model.physicsMode,'rlc')
    result.switchCurrent = x(model.index.switchCurrent,:);
    result.switchCapVoltage = x(model.index.switchCapVoltage,:);
else
    result.switchCurrent = zeros(0,nTime);
    result.switchCapVoltage = zeros(0,nTime);
end
if storeState
    result.state = x;
else
    result.state = [];
end

end

function values = local_source_values(source,t)
if isa(source,'function_handle')
    try
        values = source(t);
    catch
        values = arrayfun(source,t);
    end
elseif isnumeric(source)
    values = source;
else
    error('otr_simulate_circuit:InvalidSource', ...
        'SOURCE must be numeric or a function handle.');
end
if isscalar(values)
    values = repmat(double(values),size(t));
elseif isnumeric(values) && isreal(values) && isvector(values) && ...
        numel(values)==numel(t)
    values = double(values(:).');
else
    error('otr_simulate_circuit:InvalidSource', ...
        'SOURCE must evaluate to one real sample per time point.');
end
if any(~isfinite(values))
    error('otr_simulate_circuit:InvalidSource', ...
        'SOURCE contains non-finite values.');
end
end

function schedule = local_validate_schedule(schedule,t,nSwitch)
if ~isstruct(schedule) || ~isscalar(schedule) || ...
        ~isfield(schedule,'weights') || ~isfield(schedule,'time')
    error('otr_simulate_circuit:InvalidSchedule', ...
        'SCHEDULE must be returned by otr_switch_schedule.');
end
if numel(schedule.time)~=numel(t) || ...
        any(abs(double(schedule.time(:).')-t) > ...
        64*eps(max(1,max(abs(t)))))
    error('otr_simulate_circuit:InvalidSchedule', ...
        'The schedule and integration time grids must be identical.');
end
w = double(schedule.weights);
if size(w,2)~=numel(t)
    error('otr_simulate_circuit:InvalidSchedule', ...
        'Schedule weights must have one column per time sample.');
end
if size(w,1)==1 && nSwitch>1
    w = repmat(w,nSwitch,1);
elseif size(w,1)~=nSwitch
    error('otr_simulate_circuit:InvalidSchedule', ...
        'Schedule has %d switch rows; the model requires %d.', ...
        size(w,1),nSwitch);
end
if any(~isfinite(w(:))) || any(w(:)<0) || any(w(:)>1)
    error('otr_simulate_circuit:InvalidSchedule', ...
        'Schedule weights must lie in [0,1].');
end
schedule.weights = w;
schedule.nSwitch = nSwitch;
end

function [C,G] = local_matrices(model,weight)
C = model.Cbase;
G = model.Gbase;
if strcmp(model.physicsMode,'ideal-capacitance')
    diagonalAddition = accumarray(model.switchNodes(:), ...
        model.Cload(:).*weight(:),[model.nNodes,1]);
    node = model.index.voltage;
    C(node,node) = C(node,node) + ...
        spdiags(diagonalAddition,0,model.nNodes,model.nNodes);
else
    for k = 1:model.nSwitches
        node = model.switchNodes(k);
        branch = model.index.switchCurrent(k);
        G(node,branch) = G(node,branch)+weight(k);
        G(branch,node) = G(branch,node)-weight(k);
    end
end
end

function [xPlus,diagnostic] = local_apply_boundary(model,xMinus,oldW,newW)
xPlus = xMinus;
diagnostic = struct('relativeResidual',0,'chargeResidual',0, ...
    'voltageContinuityResidual',0,'clearedStateResidual',0, ...
    'description','no state jump');
vIndex = model.index.voltage;
vMinus = xMinus(vIndex);

if strcmp(model.physicsMode,'ideal-capacitance')
    decreasing = newW < oldW;
    increasing = newW > oldW;
    intermediateW = oldW;
    intermediateW(decreasing) = newW(decreasing);
    [Cintermediate,~] = local_matrices(model,intermediateW);
    [Cnew,~] = local_matrices(model,newW);
    Cvi = Cintermediate(vIndex,vIndex);
    Cvn = Cnew(vIndex,vIndex);
    if any(increasing)
        vPlus = Cvn \ (Cvi*vMinus);
        xPlus(vIndex) = vPlus;
        chargeResidualVector = Cvn*vPlus-Cvi*vMinus;
        scale = max([norm(Cvi*vMinus,inf),norm(Cvn*vPlus,inf),realmin]);
        chargeResidual = norm(chargeResidualVector,inf)/scale;
    else
        vPlus = vMinus;
        chargeResidual = 0;
    end
    if any(decreasing) && ~any(increasing)
        voltageResidual = norm(vPlus-vMinus,inf)/ ...
            max([norm(vMinus,inf),norm(vPlus,inf),realmin]);
    elseif any(decreasing)
        % In a mixed event, removal is applied first with no intermediate
        % voltage jump; the later insertion may change every node voltage.
        voltageResidual = 0;
    else
        voltageResidual = 0;
    end
    diagnostic.chargeResidual = chargeResidual;
    diagnostic.voltageContinuityResidual = voltageResidual;
    diagnostic.relativeResidual = max(chargeResidual,voltageResidual);
    diagnostic.description = ...
        'OFF: voltage continuous; ON: nodal conductor charge conserved';
else
    offBranches = find(newW==0 & oldW>0);
    onBranches = find(oldW==0 & newW>0);
    if ~isempty(offBranches)
        xPlus(model.index.switchCurrent(offBranches)) = 0;
        xPlus(model.index.switchCapVoltage(offBranches)) = 0;
    end
    if ~isempty(onBranches)
        % OFF branches are grounded/discharged by the physical RF switch.
        xPlus(model.index.switchCurrent(onBranches)) = 0;
        xPlus(model.index.switchCapVoltage(onBranches)) = 0;
    end
    voltageResidual = norm(xPlus(vIndex)-vMinus,inf)/ ...
        max([norm(vMinus,inf),norm(xPlus(vIndex),inf),realmin]);
    cleared = [xPlus(model.index.switchCurrent([offBranches;onBranches])); ...
        xPlus(model.index.switchCapVoltage([offBranches;onBranches]))];
    if isempty(cleared)
        clearedResidual = 0;
    else
        clearedResidual = norm(cleared,inf);
    end
    diagnostic.voltageContinuityResidual = voltageResidual;
    diagnostic.clearedStateResidual = clearedResidual;
    diagnostic.relativeResidual = max(voltageResidual,clearedResidual);
    diagnostic.description = ...
        'RLC connection changed with node voltages continuous; removed branch cleared';
end
end

function empty = local_empty_event(nState,nSwitch)
empty = struct('time',0,'index',0,'oldWeights',zeros(nSwitch,1), ...
    'newWeights',zeros(nSwitch,1),'stateMinus',zeros(nState,1), ...
    'statePlus',zeros(nState,1),'voltageJump',[], ...
    'chargeResidual',0,'voltageContinuityResidual',0, ...
    'clearedStateResidual',0,'relativeResidual',0,'description','');
end

function event = local_make_event(time,index,oldW,newW,xMinus,xPlus,d)
event = struct('time',time,'index',index,'oldWeights',oldW, ...
    'newWeights',newW,'stateMinus',xMinus,'statePlus',xPlus, ...
    'voltageJump',xPlus-xMinus,'chargeResidual',d.chargeResidual, ...
    'voltageContinuityResidual',d.voltageContinuityResidual, ...
    'clearedStateResidual',d.clearedStateResidual, ...
    'relativeResidual',d.relativeResidual,'description',d.description);
end

function value = local_field(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function value = local_logical(s,name,defaultValue)
value = local_field(s,name,defaultValue);
if ~(islogical(value) && isscalar(value)) && ...
        ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        (value==0 || value==1))
    error('otr_simulate_circuit:InvalidOption', ...
        '%s must be a logical scalar.',name);
end
value = logical(value);
end
