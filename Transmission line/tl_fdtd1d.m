function field = tl_fdtd1d(model,fdtdCfg)
%TL_FDTD1D Leapfrog simulation of a finite time-varying LC transmission line.
%
%   field = tl_fdtd1d(model,fdtdCfg)
%
% Integer-time node states are charge Q (and, for a finite Cblock, Qblock).
% Half-time branch states are current I, with Phi=Lseries*I recovered for
% output.  This preserves charge/flux continuity at capacitance time steps.
%
% Required fdtdCfg fields:
%   dt, nSteps
%
% Common optional fields:
%   boundaryType       'open' (default), 'matched', 'periodic', or 'short'
%   recordEvery        positive integer (default 1)
%   precision          'double' (default) or 'single' for histories
%   modulationEnabled  logical (default true)
%   modulationStart    integer-grid time (default 0)
%   modulationEnd      integer-grid time or Inf (default Inf)
%   nodeIndices, branchIndices
%   V0, Q0, IHalf0, PhiHalf0, I0Half0, Qblock0
%   source.type        'none', 'current', or 'thevenin'
%   source.node        node index
%   source.waveformFcn function handle returning A or V versus time
%   source.impedance   positive ohms for a Thevenin source
%   stabilityAudit     exact audit returned by an identical model/topology
%                      run; bare frequency estimates are never accepted
%
% IHalf0/PhiHalf0 are states at t=-dt/2.  IHalf is reported at t+dt/2;
% IAtNodeTime is the centered average of the adjacent half-time currents.
% The requested final state is recorded even when recordEvery does not
% divide nSteps.

if nargin ~= 2 || ~isstruct(model) || ~isscalar(model) || ...
        ~isstruct(fdtdCfg) || ~isscalar(fdtdCfg)
    error('Use tl_fdtd1d(model,fdtdCfg).');
end
required = {'kind','cell','series','shunt','resonator','ports', ...
    'finite','modulation','functions','snapshot','derived'};
missing = required(~isfield(model,required));
if ~isempty(missing)
    error('model is missing field(s): %s.',strjoin(missing,', '));
end

dt = read_positive(fdtdCfg,'dt',[]);
nSteps = read_integer(fdtdCfg,'nSteps',[],1);
recordEvery = read_integer(fdtdCfg,'recordEvery',1,1);
boundaryType = read_text(fdtdCfg,'boundaryType','open', ...
    {'open','matched','periodic','short'});
precision = read_text(fdtdCfg,'precision','double',{'double','single'});
modulationEnabled = read_logical(fdtdCfg,'modulationEnabled',true);
modulationStart = read_nonnegative(fdtdCfg,'modulationStart',0);
modulationEnd = read_positive_or_inf(fdtdCfg,'modulationEnd',Inf);
if modulationEnd <= modulationStart
    error('fdtdCfg.modulationEnd must be greater than modulationStart.');
end
assert_on_time_grid(modulationStart,dt,'modulationStart');
if isfinite(modulationEnd)
    assert_on_time_grid(modulationEnd,dt,'modulationEnd');
end

nNode = model.finite.cellCount;
isPeriodic = strcmp(boundaryType,'periodic');
if isPeriodic
    nBranch = nNode;
else
    nBranch = nNode-1;
end
[incidence,branchStart,branchEnd] = make_incidence(nNode,isPeriodic);
[Lseries,Lfactor] = make_inductance(model,nBranch,isPeriodic);

% 有限链最坏电容状态的稠密本征值审计只依赖模型和边界拓扑，不依赖
% 初值或 k_c。逐 k 波包扫描可以复用一个带完整模型快照的审计对象；
% 裸标量绝不被接受，避免以性能优化为名绕过稳定性核查。
if isfield(fdtdCfg,'stabilityAudit') && ~isempty(fdtdCfg.stabilityAudit)
    [omegaMaximum,stabilityAudit] = validate_stability_audit( ...
        fdtdCfg.stabilityAudit,model,isPeriodic,nNode,nBranch);
    stabilityAuditSource = 'validated-reuse';
else
    omegaMaximum = finite_omega_maximum(model,incidence,Lseries);
    stabilityAudit = make_stability_audit( ...
        omegaMaximum,model,isPeriodic,nNode,nBranch);
    stabilityAuditSource = 'computed-finite-chain-eigenspectrum';
end
stabilityNumber = dt*omegaMaximum;
if stabilityNumber > 1.8*(1+20*eps)
    error(['FDTD stability check failed: dt*omegaMaximum = %.6g exceeds ' ...
        'the registered limit 1.8.  Use dt <= %.6g s.'], ...
        stabilityNumber,1.8/omegaMaximum);
end

source = parse_source(fdtdCfg,model,nNode,boundaryType);
[Gport,Gsource] = terminal_conductance(model,nNode,boundaryType,source);
Gphysical = model.shunt.Gp*ones(nNode,1);
Gtotal = Gphysical+Gport+Gsource;

indices = (1:nNode).';
Cinitial = capacitance_at(0,indices,model,modulationEnabled, ...
    modulationStart,modulationEnd);
[q,currentHalf,i0Half,qblock] = initial_state( ...
    fdtdCfg,model,Cinitial,Lfactor,nNode,nBranch);
if strcmp(boundaryType,'short')
    q([1 end]) = 0;
end

Rs = model.series.Rs;
branchPlus = Lseries+(dt*Rs/2)*speye(nBranch);
branchMinus = Lseries-(dt*Rs/2)*speye(nBranch);
branchFactor = chol(branchPlus);
if strcmp(model.kind,'crow')
    L0 = model.resonator.L0;
    R0 = model.resonator.R0;
    resonatorDenominator = L0+dt*R0/2;
    resonatorNumerator = L0-dt*R0/2;
else
    resonatorDenominator = NaN;
    resonatorNumerator = NaN;
end

nodeSelection = read_indices(fdtdCfg,'nodeIndices',nNode);
branchSelection = read_indices(fdtdCfg,'branchIndices',nBranch);
recordSteps = unique([0:recordEvery:nSteps,nSteps]);
nRecord = numel(recordSteps);
history = allocate_history(nRecord,numel(nodeSelection), ...
    numel(branchSelection),model,precision);

energyTotal = zeros(1,nRecord);
energyElectric = zeros(1,nRecord);
energyMagnetic = zeros(1,nRecord);
cumulativeSource = zeros(1,nRecord);
cumulativePump = zeros(1,nRecord);
cumulativeSeriesLoss = zeros(1,nRecord);
cumulativeShuntLoss = zeros(1,nRecord);
cumulativePortLoss = zeros(1,nRecord);
sourceWork = 0;
pumpWork = 0;
seriesLoss = 0;
shuntLoss = 0;
portLoss = 0;

recordIndex = 1;
for step = 0:nSteps
    time = step*dt;
    capacitance = capacitance_at(time,indices,model,modulationEnabled, ...
        modulationStart,modulationEnd);
    voltage = q./capacitance;

    branchRhs = branchMinus*currentHalf-dt*(incidence'*voltage);
    currentPlus = solve_chol(branchFactor,branchRhs);
    currentCentered = (currentHalf+currentPlus)/2;
    phiPlus = Lseries*currentPlus;

    if strcmp(model.kind,'crow')
        if isinf(model.resonator.Cblock)
            vblock = zeros(nNode,1);
        else
            vblock = qblock/model.resonator.Cblock;
        end
        i0Plus = (resonatorNumerator*i0Half+ ...
            dt*(voltage-vblock))/resonatorDenominator;
        i0Centered = (i0Half+i0Plus)/2;
        phi0Plus = model.resonator.L0*i0Plus;
    else
        i0Plus = zeros(0,1);
        i0Centered = zeros(0,1);
        phi0Plus = zeros(0,1);
        vblock = zeros(0,1);
    end

    if recordIndex <= nRecord && step == recordSteps(recordIndex)
        history.nodeV(recordIndex,:) = cast(voltage(nodeSelection),precision);
        history.nodeQ(recordIndex,:) = cast(q(nodeSelection),precision);
        history.branchIHalf(recordIndex,:) = ...
            cast(currentPlus(branchSelection),precision);
        history.branchPhiHalf(recordIndex,:) = ...
            cast(phiPlus(branchSelection),precision);
        history.branchIAtNode(recordIndex,:) = ...
            cast(currentCentered(branchSelection),precision);
        if strcmp(model.kind,'crow')
            history.resonatorIHalf(recordIndex,:) = ...
                cast(i0Plus(nodeSelection),precision);
            history.resonatorPhiHalf(recordIndex,:) = ...
                cast(phi0Plus(nodeSelection),precision);
            if ~isinf(model.resonator.Cblock)
                history.resonatorQblock(recordIndex,:) = ...
                    cast(qblock(nodeSelection),precision);
                history.resonatorVblock(recordIndex,:) = ...
                    cast(vblock(nodeSelection),precision);
            end
        end
        [energyTotal(recordIndex),energyElectric(recordIndex), ...
            energyMagnetic(recordIndex)] = circuit_energy( ...
            q,capacitance,currentCentered,Lseries,model, ...
            i0Centered,qblock);
        cumulativeSource(recordIndex) = sourceWork;
        cumulativePump(recordIndex) = pumpWork;
        cumulativeSeriesLoss(recordIndex) = seriesLoss;
        cumulativeShuntLoss(recordIndex) = shuntLoss;
        cumulativePortLoss(recordIndex) = portLoss;
        recordIndex = recordIndex+1;
    end

    if step == nSteps
        break;
    end

    nextTime = (step+1)*dt;
    nextCapacitance = capacitance_at(nextTime,indices,model, ...
        modulationEnabled,modulationStart,modulationEnd);
    sourceCurrent = source_current(source,time+dt/2);
    forcing = zeros(nNode,1);
    forcing(source.node) = sourceCurrent;
    divergence = incidence*currentPlus;
    if strcmp(model.kind,'crow')
        divergence = divergence-i0Plus;
    end

    rhsCharge = q+dt*(divergence+forcing-0.5*Gtotal.*voltage);
    qNext = rhsCharge./(1+0.5*dt*Gtotal./nextCapacitance);
    if strcmp(boundaryType,'short')
        qNext([1 end]) = 0;
    end
    nextVoltage = qNext./nextCapacitance;
    if strcmp(model.kind,'crow') && ~isinf(model.resonator.Cblock)
        qblockNext = qblock+dt*i0Plus;
    else
        qblockNext = qblock;
    end

    voltageCentered = (voltage+nextVoltage)/2;
    sourceWork = sourceWork+dt*real(conj(voltageCentered(source.node))* ...
        sourceCurrent);
    pumpWork = pumpWork+0.5*sum(abs(q).^2.* ...
        (1./nextCapacitance-1./capacitance));
    seriesLoss = seriesLoss+dt*Rs*sum(abs(currentCentered).^2);
    shuntLoss = shuntLoss+dt*sum(Gphysical.*abs(voltageCentered).^2);
    portLoss = portLoss+dt*sum((Gport+Gsource).* ...
        abs(voltageCentered).^2);

    q = qNext;
    currentHalf = currentPlus;
    i0Half = i0Plus;
    qblock = qblockNext;
end

nodeX = (0:nNode-1)*model.cell.a;
branchX = ((branchStart+branchEnd)/2-1)*model.cell.a;
if isPeriodic
    branchX(end) = (nNode-0.5)*model.cell.a;
end
nodeTime = recordSteps*dt;

field.node = struct('x',nodeX(nodeSelection), ...
    't',nodeTime,'V',history.nodeV,'Q',history.nodeQ, ...
    'indices',nodeSelection);
field.branch = struct('x',branchX(branchSelection), ...
    'tHalf',nodeTime+dt/2,'IHalf',history.branchIHalf, ...
    'PhiHalf',history.branchPhiHalf, ...
    'IAtNodeTime',history.branchIAtNode, ...
    'indices',branchSelection,'startNode',branchStart(branchSelection), ...
    'endNode',branchEnd(branchSelection));
if strcmp(model.kind,'crow')
    field.resonator = struct('x',nodeX(nodeSelection),'tHalf',nodeTime+dt/2, ...
        'IHalf',history.resonatorIHalf, ...
        'PhiHalf',history.resonatorPhiHalf);
    if ~isinf(model.resonator.Cblock)
        field.resonator.Qblock = history.resonatorQblock;
        field.resonator.Vblock = history.resonatorVblock;
        field.resonator.t = nodeTime;
    end
end

energyScale = max([abs(energyTotal),abs(cumulativeSource), ...
    abs(cumulativePump),abs(cumulativeSeriesLoss), ...
    abs(cumulativeShuntLoss),abs(cumulativePortLoss),realmin]);
ledgerResidual = (energyTotal-energyTotal(1))-cumulativeSource- ...
    cumulativePump+cumulativeSeriesLoss+cumulativeShuntLoss+ ...
    cumulativePortLoss;
field.energy = struct('t',nodeTime,'total',energyTotal, ...
    'electric',energyElectric,'magnetic',energyMagnetic, ...
    'sourceWork',cumulativeSource,'pumpWork',cumulativePump, ...
    'seriesDissipation',cumulativeSeriesLoss, ...
    'shuntDissipation',cumulativeShuntLoss, ...
    'portDissipation',cumulativePortLoss,'ledgerResidual',ledgerResidual, ...
    'relativeLedgerResidual',ledgerResidual/energyScale);
field.grid = struct('nodeCount',nNode,'branchCount',nBranch, ...
    'dt',dt,'nSteps',nSteps,'recordEvery',recordEvery, ...
    'recordSteps',recordSteps,'boundaryType',boundaryType, ...
    'stabilityNumber',stabilityNumber,'omegaMaximum',omegaMaximum, ...
    'stabilityAuditSource',stabilityAuditSource, ...
    'stabilityAudit',stabilityAudit);
field.source = source_snapshot(source);
field.modulation = struct('enabled',modulationEnabled, ...
    'start',modulationStart,'end',modulationEnd);
field.modelSnapshot = model.snapshot;
field.fdtdConfig = sanitize_config(fdtdCfg);
field.observableNote = ['I is a transmission-line proxy for local magnetic ' ...
    'field; V is a proxy for local electric field, not a PCB full-wave field.'];
end

% -------------------------------------------------------------------------
function [B,startNode,endNode] = make_incidence(nNode,isPeriodic)
if isPeriodic
    startNode = 1:nNode;
    endNode = [2:nNode,1];
else
    startNode = 1:nNode-1;
    endNode = 2:nNode;
end
nBranch = numel(startNode);
rows = [startNode,endNode];
columns = [1:nBranch,1:nBranch];
values = [-ones(1,nBranch),ones(1,nBranch)];
B = sparse(rows,columns,values,nNode,nBranch);
end

% -------------------------------------------------------------------------
function [Lseries,Lfactor] = make_inductance(model,nBranch,isPeriodic)
Ls = model.series.Ls;
S = model.series.mutualS;
Lseries = spdiags(Ls*ones(nBranch,1),0,nBranch,nBranch);
if nBranch > 1 && S ~= 0
    Lseries = Lseries+spdiags(S*ones(nBranch,1),1,nBranch,nBranch)+ ...
        spdiags(S*ones(nBranch,1),-1,nBranch,nBranch);
    if isPeriodic && nBranch > 2
        Lseries(1,end) = S;
        Lseries(end,1) = S;
    end
end
[Lfactor,flag] = chol(Lseries);
if flag ~= 0
    error('The finite-chain mutual-inductance matrix is not positive definite.');
end
end

% -------------------------------------------------------------------------
function omegaMaximum = finite_omega_maximum(model,B,Lseries)
nNode = size(B,1);
nBranch = size(B,2);
scale = model.finite.amplitudeScaleByCell;
if isscalar(scale), scale = repmat(scale,nNode,1); end
Cminimum = model.shunt.totalC0-model.shunt.deltaC*scale(:);
Cinv = spdiags(1./Cminimum,0,nNode,nNode);
LinvBt = Lseries\B';
if strcmp(model.kind,'sspp')
    A = [sparse(nNode,nNode),B; -LinvBt*Cinv,sparse(nBranch,nBranch)];
else
    n0 = nNode;
    top = [sparse(nNode,nNode),B,-speye(nNode)];
    series = [-LinvBt*Cinv,sparse(nBranch,nBranch), ...
        sparse(nBranch,n0)];
    if isinf(model.resonator.Cblock)
        resonator = [Cinv/model.resonator.L0, ...
            sparse(n0,nBranch),sparse(n0,n0)];
        A = [top;series;resonator];
    else
        top = [top,sparse(nNode,n0)];
        series = [series,sparse(nBranch,n0)];
        resonator = [Cinv/model.resonator.L0, ...
            sparse(n0,nBranch),sparse(n0,n0), ...
            -speye(n0)/(model.resonator.L0*model.resonator.Cblock)];
        blockCharge = [sparse(n0,nNode+nBranch),speye(n0),sparse(n0,n0)];
        A = [top;series;resonator;blockCharge];
    end
end
roots = eig(full(A));
omegaMaximum = max(abs(roots));
if ~isfinite(omegaMaximum) || omegaMaximum <= 0
    error('Unable to determine a positive finite-chain frequency bound.');
end
end

% -------------------------------------------------------------------------
function audit = make_stability_audit( ...
        omegaMaximum,model,isPeriodic,nNode,nBranch)
audit = struct();
audit.version = 'tl-fdtd-stability-v1';
audit.modelSnapshot = model.snapshot;
audit.isPeriodic = isPeriodic;
audit.nodeCount = nNode;
audit.branchCount = nBranch;
audit.omegaMaximum = omegaMaximum;
audit.method = ['finite-chain conservative eigenfrequency bound at the ' ...
    'registered per-node minimum capacitances'];
end

% -------------------------------------------------------------------------
function [omegaMaximum,audit] = validate_stability_audit( ...
        value,model,isPeriodic,nNode,nBranch)
required = {'version','modelSnapshot','isPeriodic','nodeCount', ...
    'branchCount','omegaMaximum','method'};
if ~isstruct(value) || ~isscalar(value) || ...
        any(~isfield(value,required)) || ...
        ~strcmp(value.version,'tl-fdtd-stability-v1') || ...
        ~islogical(value.isPeriodic) || ~isscalar(value.isPeriodic) || ...
        value.isPeriodic ~= isPeriodic || ...
        ~isequal(value.nodeCount,nNode) || ...
        ~isequal(value.branchCount,nBranch) || ...
        ~isequaln(value.modelSnapshot,model.snapshot) || ...
        ~isnumeric(value.omegaMaximum) || ...
        ~isscalar(value.omegaMaximum) || ...
        ~isreal(value.omegaMaximum) || ...
        ~isfinite(value.omegaMaximum) || value.omegaMaximum <= 0
    error(['fdtdCfg.stabilityAudit does not exactly match the current ' ...
        'model and boundary topology; omit it to recompute the audit.']);
end
omegaMaximum = value.omegaMaximum;
audit = value;
end

% -------------------------------------------------------------------------
function source = parse_source(cfg,model,nNode,boundaryType)
if ~isfield(cfg,'source') || isempty(cfg.source)
    input = struct('type','none');
else
    input = cfg.source;
end
if ~isstruct(input) || ~isscalar(input)
    error('fdtdCfg.source must be a scalar struct.');
end
source.type = read_text(input,'type','none',{'none','current','thevenin'});
source.node = read_integer(input,'node',1,1);
if source.node > nNode
    error('fdtdCfg.source.node exceeds the finite node count.');
end
if strcmp(boundaryType,'short') && any(source.node == [1 nNode]) && ...
        ~strcmp(source.type,'none')
    error('A source cannot drive a node constrained by a short boundary.');
end
if strcmp(source.type,'none')
    source.waveformFcn = @(t) 0*t;
    source.impedance = Inf;
else
    if ~isfield(input,'waveformFcn') || ...
            ~isa(input.waveformFcn,'function_handle')
        error('An active source requires source.waveformFcn.');
    end
    source.waveformFcn = input.waveformFcn;
    if strcmp(source.type,'thevenin')
        source.impedance = read_positive(input,'impedance',model.ports.Zsource);
    else
        source.impedance = Inf;
    end
end
end

% -------------------------------------------------------------------------
function [Gport,Gsource] = terminal_conductance(model,nNode,boundary,source)
Gport = zeros(nNode,1);
Gsource = zeros(nNode,1);
if strcmp(boundary,'matched')
    Gport(1) = 1/model.ports.Zsource;
    Gport(end) = Gport(end)+1/model.ports.Zload;
end
if strcmp(source.type,'thevenin')
    alreadyPresent = strcmp(boundary,'matched') && source.node == 1 && ...
        abs(source.impedance-model.ports.Zsource) <= ...
        100*eps(max(source.impedance,model.ports.Zsource));
    if ~alreadyPresent
        Gsource(source.node) = 1/source.impedance;
    end
end
end

% -------------------------------------------------------------------------
function value = source_current(source,time)
sample = source.waveformFcn(time);
if ~isnumeric(sample) || ~isscalar(sample) || ~isfinite(sample)
    error('source.waveformFcn must return one finite numeric scalar.');
end
if strcmp(source.type,'thevenin')
    value = sample/source.impedance;
elseif strcmp(source.type,'current')
    value = sample;
else
    value = 0;
end
end

% -------------------------------------------------------------------------
function C = capacitance_at(time,indices,model,enabled,startTime,endTime)
active = enabled && time >= startTime && time < endTime;
if active
    C = model.functions.capacitanceFinite(indices,time-startTime);
else
    C = model.shunt.totalC0*ones(numel(indices),1);
end
C = C(:);
end

% -------------------------------------------------------------------------
function [q,currentHalf,i0Half,qblock] = initial_state( ...
        cfg,model,Cinitial,Lfactor,nNode,nBranch)
hasQ = isfield(cfg,'Q0') && ~isempty(cfg.Q0);
hasV = isfield(cfg,'V0') && ~isempty(cfg.V0);
if hasQ && hasV
    error('Specify only one of fdtdCfg.Q0 and fdtdCfg.V0.');
elseif hasQ
    q = read_state(cfg.Q0,nNode,'Q0');
elseif hasV
    q = Cinitial.*read_state(cfg.V0,nNode,'V0');
else
    q = zeros(nNode,1);
end
hasI = isfield(cfg,'IHalf0') && ~isempty(cfg.IHalf0);
hasPhi = isfield(cfg,'PhiHalf0') && ~isempty(cfg.PhiHalf0);
if hasI && hasPhi
    error('Specify only one of fdtdCfg.IHalf0 and fdtdCfg.PhiHalf0.');
elseif hasI
    currentHalf = read_state(cfg.IHalf0,nBranch,'IHalf0');
elseif hasPhi
    currentHalf = solve_chol(Lfactor,read_state( ...
        cfg.PhiHalf0,nBranch,'PhiHalf0'));
else
    currentHalf = zeros(nBranch,1);
end
if strcmp(model.kind,'crow')
    if isfield(cfg,'I0Half0') && ~isempty(cfg.I0Half0)
        i0Half = read_state(cfg.I0Half0,nNode,'I0Half0');
    else
        i0Half = zeros(nNode,1);
    end
    if ~isinf(model.resonator.Cblock)
        if isfield(cfg,'Qblock0') && ~isempty(cfg.Qblock0)
            qblock = read_state(cfg.Qblock0,nNode,'Qblock0');
        else
            qblock = zeros(nNode,1);
        end
    else
        qblock = zeros(0,1);
    end
else
    i0Half = zeros(0,1);
    qblock = zeros(0,1);
end
end

% -------------------------------------------------------------------------
function state = read_state(value,count,label)
if ~isnumeric(value) || ~isvector(value) || any(~isfinite(value)) || ...
        ~(isscalar(value) || numel(value) == count)
    error('fdtdCfg.%s must be a finite scalar or have the required size.',label);
end
if isscalar(value)
    state = repmat(value,count,1);
else
    state = value(:);
end
end

% -------------------------------------------------------------------------
function selection = read_indices(cfg,name,count)
if ~isfield(cfg,name) || isempty(cfg.(name))
    selection = 1:count;
else
    selection = cfg.(name);
    if ~isnumeric(selection) || ~isvector(selection) || ...
            any(~isfinite(selection)) || any(selection ~= round(selection)) || ...
            any(selection < 1) || any(selection > count) || ...
            numel(unique(selection)) ~= numel(selection)
        error('fdtdCfg.%s must contain unique valid indices.',name);
    end
    selection = selection(:).';
end
end

% -------------------------------------------------------------------------
function history = allocate_history(nRecord,nNode,nBranch,model,precision)
history.nodeV = zeros(nRecord,nNode,precision);
history.nodeQ = zeros(nRecord,nNode,precision);
history.branchIHalf = zeros(nRecord,nBranch,precision);
history.branchPhiHalf = zeros(nRecord,nBranch,precision);
history.branchIAtNode = zeros(nRecord,nBranch,precision);
if strcmp(model.kind,'crow')
    history.resonatorIHalf = zeros(nRecord,nNode,precision);
    history.resonatorPhiHalf = zeros(nRecord,nNode,precision);
    if ~isinf(model.resonator.Cblock)
        history.resonatorQblock = zeros(nRecord,nNode,precision);
        history.resonatorVblock = zeros(nRecord,nNode,precision);
    end
end
end

% -------------------------------------------------------------------------
function [total,electric,magnetic] = circuit_energy( ...
        q,C,current,Lseries,model,i0,qblock)
electric = 0.5*sum(abs(q).^2./C);
magnetic = 0.5*real(current'*(Lseries*current));
if strcmp(model.kind,'crow')
    magnetic = magnetic+0.5*model.resonator.L0*sum(abs(i0).^2);
    if ~isinf(model.resonator.Cblock)
        electric = electric+0.5*sum(abs(qblock).^2)/model.resonator.Cblock;
    end
end
total = electric+magnetic;
end

% -------------------------------------------------------------------------
function solution = solve_chol(factor,rhs)
solution = factor\(factor'\rhs);
end

% -------------------------------------------------------------------------
function snapshot = source_snapshot(source)
snapshot = source;
snapshot = rmfield(snapshot,'waveformFcn');
snapshot.waveformDescription = 'function handle omitted from snapshot';
end

% -------------------------------------------------------------------------
function snapshot = sanitize_config(cfg)
snapshot = cfg;
if isfield(snapshot,'source') && isstruct(snapshot.source) && ...
        isfield(snapshot.source,'waveformFcn')
    snapshot.source = rmfield(snapshot.source,'waveformFcn');
    snapshot.source.waveformDescription = 'function handle omitted from snapshot';
end
end

% -------------------------------------------------------------------------
function assert_on_time_grid(value,dt,label)
gridIndex = value/dt;
if abs(gridIndex-round(gridIndex)) > 100*eps(max(1,abs(gridIndex)))
    error('fdtdCfg.%s must lie on an integer FDTD time node.',label);
end
end

% -------------------------------------------------------------------------
function value = read_text(cfg,name,defaultValue,allowed)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || size(value,1) ~= 1
    error('%s must be a text scalar.',name);
end
match = strcmpi(value,allowed);
if ~any(match), error('%s must be one of: %s.',name,strjoin(allowed,', ')); end
value = allowed{find(match,1)};
end

function value = read_logical(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~islogical(value) || ~isscalar(value), error('%s must be a logical scalar.',name); end
end

function value = read_positive(cfg,name,defaultValue)
value = read_real(cfg,name,defaultValue);
if value <= 0, error('%s must be positive.',name); end
end

function value = read_nonnegative(cfg,name,defaultValue)
value = read_real(cfg,name,defaultValue);
if value < 0, error('%s must be nonnegative.',name); end
end

function value = read_positive_or_inf(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || isnan(value) || value <= 0
    error('%s must be positive and may be Inf.',name);
end
end

function value = read_integer(cfg,name,defaultValue,minimum)
value = read_real(cfg,name,defaultValue);
if value ~= round(value) || value < minimum
    error('%s must be an integer not smaller than %d.',name,minimum);
end
end

function value = read_real(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(value)
    error('%s must be a finite real scalar.',name);
end
end
