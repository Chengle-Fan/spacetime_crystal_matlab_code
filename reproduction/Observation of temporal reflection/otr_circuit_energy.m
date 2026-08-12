function energy = otr_circuit_energy(model, result)
%OTR_CIRCUIT_ENERGY Audit storage, loss, source work, and switch work.
%   ENERGY = OTR_CIRCUIT_ENERGY(MODEL,RESULT) evaluates the energy balance
%   of a result returned by OTR_SIMULATE_CIRCUIT.  All quantities are in SI
%   units.  The stored energy is formed directly from physical capacitors
%   and inductors; dissipation includes the line/shunt losses, source
%   resistance, 50-ohm load, and RLC switch losses where applicable.
%
%   Instantaneous switch events are included as an explicit jump work
%
%       W_switch(t_e) = E(t_e^+) - E(t_e^-).
%
%   A negative value includes charge-sharing loss or energy carried away by
%   a removed capacitor.  Between switch events, trapezoidal MNA has the
%   exact discrete quadratic balance when power is evaluated at its
%   arithmetic midpoint.  The reported residual is therefore also a useful
%   solver and boundary-condition test.

if ~isstruct(model) || ~isscalar(model) || ~isfield(model,'isMNA') || ...
        ~model.isMNA
    error('otr_circuit_energy:InvalidModel', ...
        'MODEL must be returned by otr_build_mna.');
end
required = {'time','source','weights','voltage','lineCurrent','eventRecords'};
if ~isstruct(result) || ~isscalar(result) || ...
        ~all(isfield(result,required))
    error('otr_circuit_energy:InvalidResult', ...
        'RESULT must be returned by otr_simulate_circuit.');
end
t = double(result.time(:).');
nTime = numel(t);
if nTime<2 || any(diff(t)<=0)
    error('otr_circuit_energy:InvalidResult', ...
        'RESULT contains an invalid time grid.');
end
if size(result.voltage,2)~=nTime || size(result.voltage,1)~=model.nNodes || ...
        size(result.lineCurrent,2)~=nTime || ...
        size(result.weights,2)~=nTime
    error('otr_circuit_energy:InvalidResult', ...
        'RESULT state arrays have inconsistent sizes.');
end

state = local_reconstruct_state(model,result);
stored = zeros(1,nTime);
electric = zeros(1,nTime);
magnetic = zeros(1,nTime);
for sample = 1:nTime
    [stored(sample),electric(sample),magnetic(sample)] = ...
        local_state_energy(model,state(:,sample),result.weights(:,sample));
end

% The right-hand state used during each integration interval is x(t_n^-),
% not x(t_n^+) when a switch event exists at t_n.
rightMinus = state;
eventJump = zeros(1,nTime);
eventTable = result.eventRecords;
for event = 1:numel(eventTable)
    index = eventTable(event).index;
    if index<1 || index>nTime
        error('otr_circuit_energy:InvalidResult', ...
            'An event index lies outside the result time grid.');
    end
    rightMinus(:,index) = eventTable(event).stateMinus;
    [minusEnergy,~,~] = local_state_energy(model, ...
        eventTable(event).stateMinus,eventTable(event).oldWeights);
    [plusEnergy,~,~] = local_state_energy(model, ...
        eventTable(event).statePlus,eventTable(event).newWeights);
    eventJump(index) = eventJump(index) + plusEnergy-minusEnergy;
end

nInterval = nTime-1;
sourcePowerMid = zeros(1,nInterval);
dissipationPowerMid = zeros(1,nInterval);
loadPowerMid = zeros(1,nInterval);
sourceResistancePowerMid = zeros(1,nInterval);
lineLossPowerMid = zeros(1,nInterval);
shuntLossPowerMid = zeros(1,nInterval);
switchLossPowerMid = zeros(1,nInterval);
intervalNetWork = zeros(1,nInterval);
for interval = 1:nInterval
    h = t(interval+1)-t(interval);
    xMid = (state(:,interval)+rightMinus(:,interval+1))/2;
    uMid = (result.source(interval)+result.source(interval+1))/2;
    v = xMid(model.index.voltage);
    iLine = xMid(model.index.lineCurrent);
    inputCurrent = (uMid-v(1))/model.parameters.sourceResistance;
    sourcePowerMid(interval) = uMid*inputCurrent;
    sourceResistancePowerMid(interval) = ...
        (uMid-v(1))^2/model.parameters.sourceResistance;
    loadPowerMid(interval) = v(end)^2/model.parameters.loadResistance;
    lineLossPowerMid(interval) = sum(model.parameters.Rseries.*iLine.^2);
    shuntLossPowerMid(interval) = v.'*model.GnodePhysical*v;
    if strcmp(model.physicsMode,'rlc')
        iSwitch = xMid(model.index.switchCurrent);
        switchLossPowerMid(interval) = sum( ...
            (model.parameters.Ron+model.parameters.Rpar).*iSwitch.^2);
    end
    dissipationPowerMid(interval) = sourceResistancePowerMid(interval)+ ...
        loadPowerMid(interval)+lineLossPowerMid(interval)+ ...
        shuntLossPowerMid(interval)+switchLossPowerMid(interval);
    intervalNetWork(interval) = h*(sourcePowerMid(interval)- ...
        dissipationPowerMid(interval));
end

cumulativeContinuousWork = [0,cumsum(intervalNetWork)];
cumulativeSwitchWork = cumsum(eventJump);
balanceResidual = stored-stored(1)-cumulativeContinuousWork- ...
    cumulativeSwitchWork;
scale = max([max(abs(stored)),max(abs(cumulativeContinuousWork)), ...
    max(abs(cumulativeSwitchWork)),realmin]);

energy = struct();
energy.time = t;
energy.stored = stored;
energy.electric = electric;
energy.magnetic = magnetic;
energy.sourcePowerMid = sourcePowerMid;
energy.dissipationPowerMid = dissipationPowerMid;
energy.loadPowerMid = loadPowerMid;
energy.sourceResistancePowerMid = sourceResistancePowerMid;
energy.lineLossPowerMid = lineLossPowerMid;
energy.shuntLossPowerMid = shuntLossPowerMid;
energy.switchLossPowerMid = switchLossPowerMid;
energy.intervalContinuousWork = intervalNetWork;
energy.cumulativeContinuousWork = cumulativeContinuousWork;
energy.eventEnergyJump = eventJump;
energy.cumulativeSwitchWork = cumulativeSwitchWork;
energy.balanceResidual = balanceResidual;
energy.maxAbsBalanceResidual = max(abs(balanceResidual));
energy.relativeBalanceResidual = energy.maxAbsBalanceResidual/scale;
energy.balanceEquation = ...
    'E-E0 = integral(Psource-Pdissipation)dt + sum(Eplus-Eminus)';

end

function state = local_reconstruct_state(model,result)
if isfield(result,'state') && ~isempty(result.state)
    state = double(result.state);
    if ~isequal(size(state),[model.nState,numel(result.time)])
        error('otr_circuit_energy:InvalidResult', ...
            'RESULT.state has an inconsistent size.');
    end
    return;
end
state = zeros(model.nState,numel(result.time));
state(model.index.voltage,:) = result.voltage;
state(model.index.lineCurrent,:) = result.lineCurrent;
if strcmp(model.physicsMode,'rlc')
    if ~isfield(result,'switchCurrent') || ...
            ~isfield(result,'switchCapVoltage')
        error('otr_circuit_energy:MissingState', ...
            'RLC energy requires switchCurrent and switchCapVoltage.');
    end
    state(model.index.switchCurrent,:) = result.switchCurrent;
    state(model.index.switchCapVoltage,:) = result.switchCapVoltage;
end
end

function [total,electric,magnetic] = local_state_energy(model,x,weight)
v = x(model.index.voltage);
Cnode = model.CnodeBase;
if strcmp(model.physicsMode,'ideal-capacitance')
    addition = accumarray(model.switchNodes(:), ...
        model.Cload(:).*weight(:),[model.nNodes,1]);
    Cnode = Cnode+spdiags(addition,0,model.nNodes,model.nNodes);
end
electric = 0.5*real(v.'*Cnode*v);
iLine = x(model.index.lineCurrent);
magnetic = 0.5*sum(model.parameters.Lseries.*iLine.^2);
if strcmp(model.physicsMode,'rlc')
    iSwitch = x(model.index.switchCurrent);
    vCap = x(model.index.switchCapVoltage);
    magnetic = magnetic+0.5*sum(model.parameters.Lpar.*iSwitch.^2);
    electric = electric+0.5*sum(model.parameters.Cload.*vCap.^2);
end
total = electric+magnetic;
end
