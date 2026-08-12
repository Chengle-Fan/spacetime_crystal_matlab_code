function model = otr_build_mna(options)
%OTR_BUILD_MNA Build a Kirchhoff/MNA model of the switched transmission line.
%   MODEL = OTR_BUILD_MNA(OPTIONS) constructs a lumped, discrete circuit.
%   It is not a continuous-medium FDTD discretisation.  The unknowns are
%   node voltages and physical branch currents/voltages, and every matrix
%   row is either a KCL or a branch constitutive equation.
%
%   The default ladder contains 30 series-RL cells, 31 voltage nodes, a
%   Thevenin source at the first node, a 50 ohm load at the last node, and
%   30 switchable shunt loads at nodes 2:31.  The cell L and C are obtained
%   from Z0, d, and the microstrip effective dielectric constant.
%
%   OPTIONS.physicsMode is one of
%       'ideal-capacitance'  instantaneous capacitance insertion/removal;
%                            closing obeys C_old*v- = C_new*v+ and opening
%                            obeys v+ = v- (SI equations S14 and S15).
%       'rlc'                a fixed maximum topology with a series
%                            (Ron+Rpar)-Lpar-Cload branch at every switch
%                            node.  A connection weight couples/decouples
%                            each branch without changing the state size.
%
%   Important fields accepted in OPTIONS (all SI units) are:
%       nCells, d, Z0, epsEff, c0, Lseries, Ccell, Rseries,
%       Gshunt, sourceResistance, loadResistance, Cload, Ron, Rpar,
%       Lpar, includeParasitics, Cnearest, CnextNearest, switchNodes.
%
%   The descriptor equation returned by this routine has the form
%       C(w) * xdot + G(w) * x = B * vSource(t).
%   Use OTR_SIMULATE_CIRCUIT to assemble C(w), G(w), apply temporal
%   boundary conditions, and integrate this equation.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('otr_build_mna:InvalidOptions', ...
        'OPTIONS must be a scalar structure.');
end

% Parameters reported in the paper/Methods.
p.nCells = local_field(options, 'nCells', 30);
p.d = local_field(options, 'd', 0.2080);
p.Z0 = local_field(options, 'Z0', 50);
p.epsEff = local_field(options, 'epsEff', 8.36);
p.c0 = local_field(options, 'c0', 299792458);
p.physicsMode = lower(char(local_field(options, 'physicsMode', ...
    'ideal-capacitance')));

local_positive_integer(p.nCells, 'nCells');
local_positive_scalar(p.d, 'd');
local_positive_scalar(p.Z0, 'Z0');
local_positive_scalar(p.epsEff, 'epsEff');
local_positive_scalar(p.c0, 'c0');
if ~ismember(p.physicsMode, {'ideal-capacitance', 'rlc'})
    error('otr_build_mna:InvalidPhysicsMode', ...
        'physicsMode must be ''ideal-capacitance'' or ''rlc''.');
end

p.phaseVelocity = p.c0 / sqrt(p.epsEff);
p.cellDelay = p.d / p.phaseVelocity;
p.Lseries = local_field(options, 'Lseries', p.Z0 * p.cellDelay);
p.Ccell = local_field(options, 'Ccell', p.cellDelay / p.Z0);
p.Rseries = local_field(options, 'Rseries', 0.25);
p.lossTangent = local_field(options, 'lossTangent', 0.0019);
p.referenceFrequency = local_field(options, 'referenceFrequency', 100e6);
defaultG = 2*pi*p.referenceFrequency*p.Ccell*p.lossTangent;
p.Gshunt = local_field(options, 'Gshunt', defaultG);
p.sourceResistance = local_field(options, 'sourceResistance', 50);
p.loadResistance = local_field(options, 'loadResistance', 50);
p.Cload = local_field(options, 'Cload', 82e-12);
p.Ron = local_field(options, 'Ron', 1);
p.Rpar = local_field(options, 'Rpar', 4);
p.Lpar = local_field(options, 'Lpar', 8e-9);
p.includeParasitics = local_logical_field(options, ...
    'includeParasitics', true);
p.Cnearest = local_field(options, 'Cnearest', 2e-12);
p.CnextNearest = local_field(options, 'CnextNearest', 1e-12);

nV = p.nCells + 1;
nLine = p.nCells;
defaultSwitchNodes = 2:nV;
p.switchNodes = local_field(options, 'switchNodes', defaultSwitchNodes);
p.switchNodes = double(p.switchNodes(:).');
if isempty(p.switchNodes) || any(~isfinite(p.switchNodes)) || ...
        any(p.switchNodes ~= round(p.switchNodes)) || ...
        any(p.switchNodes < 1) || any(p.switchNodes > nV) || ...
        numel(unique(p.switchNodes)) ~= numel(p.switchNodes)
    error('otr_build_mna:InvalidSwitchNodes', ...
        'switchNodes must contain unique integer node indices in 1:%d.', nV);
end
nSwitch = numel(p.switchNodes);

p.Lseries = local_expand(p.Lseries, nLine, 'Lseries', true);
p.Rseries = local_expand(p.Rseries, nLine, 'Rseries', false);
p.Cload = local_expand(p.Cload, nSwitch, 'Cload', true);
p.Ron = local_expand(p.Ron, nSwitch, 'Ron', false);
p.Rpar = local_expand(p.Rpar, nSwitch, 'Rpar', false);
p.Lpar = local_expand(p.Lpar, nSwitch, 'Lpar', true);
local_nonnegative_scalar(p.Gshunt, 'Gshunt');
local_positive_scalar(p.sourceResistance, 'sourceResistance');
local_positive_scalar(p.loadResistance, 'loadResistance');
local_nonnegative_scalar(p.Cnearest, 'Cnearest');
local_nonnegative_scalar(p.CnextNearest, 'CnextNearest');

% Series branch incidence: current is positive from node j to node j+1.
A = sparse([1:nLine, 2:nV], [1:nLine, 1:nLine], ...
    [ones(1,nLine), -ones(1,nLine)], nV, nLine);

% A cell's distributed shunt capacitance/conductance is split equally
% between its two endpoints.  Interior nodes therefore receive one Ccell.
cNodeDiag = zeros(nV,1);
gNodeDiag = zeros(nV,1);
for branch = 1:nLine
    cNodeDiag(branch:branch+1) = cNodeDiag(branch:branch+1) + p.Ccell/2;
    gNodeDiag(branch:branch+1) = gNodeDiag(branch:branch+1) + p.Gshunt/2;
end
CnodeBase = spdiags(cNodeDiag, 0, nV, nV);

% Measured parasitic electric coupling between switch cells.  These are
% genuine capacitors between nodes, so they add the Laplacian stamp
% C*[1 -1; -1 1] to the nodal capacitance matrix.
if p.includeParasitics
    for separation = 1:2
        if separation == 1
            couplingC = p.Cnearest;
        else
            couplingC = p.CnextNearest;
        end
        for j = 1:(nSwitch-separation)
            CnodeBase = local_stamp_mutual_cap(CnodeBase, ...
                p.switchNodes(j), p.switchNodes(j+separation), couplingC);
        end
    end
end

% Resistive stamps at voltage nodes.  The Thevenin source is represented
% by its Norton equivalent: vSource/Rs in parallel with Rs.
GnodePhysical = spdiags(gNodeDiag, 0, nV, nV);
Gnode = GnodePhysical;
Gnode(1,1) = Gnode(1,1) + 1/p.sourceResistance;
Gnode(end,end) = Gnode(end,end) + 1/p.loadResistance;

nBase = nV + nLine;
Cbase = blkdiag(CnodeBase, spdiags(p.Lseries,0,nLine,nLine));
Gbase = [Gnode, A; -A.', spdiags(p.Rseries,0,nLine,nLine)];
Bbase = sparse(1,1,1/p.sourceResistance,nBase,1);

index.voltage = 1:nV;
index.lineCurrent = nV + (1:nLine);
index.switchCurrent = [];
index.switchCapVoltage = [];

switch p.physicsMode
    case 'ideal-capacitance'
        Cdescriptor = Cbase;
        Gdescriptor = Gbase;
        B = Bbase;
        nState = nBase;

    case 'rlc'
        % State ordering: [node voltages; line currents; switch branch
        % currents; switch capacitor voltages].  The fixed internal RLC
        % equations are
        %   Lp di/dt + (Ron+Rpar)i + vc - w*vnode = 0,
        %   Cload dvc/dt - i = 0.
        % KCL receives +w*i.  The +/-w coupling is skew-symmetric and is
        % therefore a lossless connection for fixed w.
        index.switchCurrent = nBase + (1:nSwitch);
        index.switchCapVoltage = nBase + nSwitch + (1:nSwitch);
        nState = nBase + 2*nSwitch;
        Cdescriptor = blkdiag(Cbase, ...
            spdiags(p.Lpar,0,nSwitch,nSwitch), ...
            spdiags(p.Cload,0,nSwitch,nSwitch));
        Gdescriptor = sparse(nState,nState);
        Gdescriptor(1:nBase,1:nBase) = Gbase;
        Gdescriptor(index.switchCurrent,index.switchCurrent) = ...
            spdiags(p.Ron+p.Rpar,0,nSwitch,nSwitch);
        Gdescriptor(index.switchCurrent,index.switchCapVoltage) = ...
            speye(nSwitch);
        Gdescriptor(index.switchCapVoltage,index.switchCurrent) = ...
            -speye(nSwitch);
        B = sparse(nState,1);
        B(1:nBase) = Bbase;
end

model = struct();
model.name = 'Discrete switched transmission-line MNA model';
model.form = 'C(w)*xdot + G(w)*x = B*vSource(t)';
model.physicsMode = p.physicsMode;
model.parameters = p;
model.nCells = p.nCells;
model.nNodes = nV;
model.nLineBranches = nLine;
model.nSwitches = nSwitch;
model.nState = nState;
model.index = index;
model.incidence = A;
model.CnodeBase = CnodeBase;
model.GnodePhysical = GnodePhysical;
model.Cbase = Cdescriptor;
model.Gbase = Gdescriptor;
model.B = B;
model.switchNodes = p.switchNodes;
model.Cload = p.Cload;
model.isMNA = true;
model.notes = { ...
    'Every voltage row is KCL; every current row is an inductor/RLC law.', ...
    'The source is a Thevenin source represented by its exact Norton stamp.', ...
    'RLC mode uses a lossless connection weight and clears a removed branch.', ...
    'Ideal-capacitance mode implements SI boundary conditions S14/S15.'};

end

function value = local_field(s, name, defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function value = local_logical_field(s, name, defaultValue)
value = local_field(s,name,defaultValue);
if ~(islogical(value) && isscalar(value)) && ...
        ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        (value == 0 || value == 1))
    error('otr_build_mna:InvalidLogical', '%s must be a logical scalar.', name);
end
value = logical(value);
end

function local_positive_integer(value, name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value >= 1 && value == round(value))
    error('otr_build_mna:InvalidParameter', ...
        '%s must be a positive integer.', name);
end
end

function local_positive_scalar(value, name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value > 0)
    error('otr_build_mna:InvalidParameter', ...
        '%s must be a positive finite scalar.', name);
end
end

function local_nonnegative_scalar(value, name)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value >= 0)
    error('otr_build_mna:InvalidParameter', ...
        '%s must be a nonnegative finite scalar.', name);
end
end

function vector = local_expand(value, n, name, strictlyPositive)
if ~(isnumeric(value) && isreal(value) && all(isfinite(value(:))))
    error('otr_build_mna:InvalidParameter', '%s must be real and finite.', name);
end
if isscalar(value)
    vector = repmat(double(value),n,1);
elseif numel(value) == n
    vector = double(value(:));
else
    error('otr_build_mna:InvalidParameter', ...
        '%s must be scalar or have %d elements.', name, n);
end
if strictlyPositive
    valid = all(vector > 0);
else
    valid = all(vector >= 0);
end
if ~valid
    qualifier = 'nonnegative';
    if strictlyPositive
        qualifier = 'positive';
    end
    error('otr_build_mna:InvalidParameter', ...
        '%s entries must be %s.', name, qualifier);
end
end

function C = local_stamp_mutual_cap(C, node1, node2, capacitance)
if capacitance == 0
    return;
end
C(node1,node1) = C(node1,node1) + capacitance;
C(node2,node2) = C(node2,node2) + capacitance;
C(node1,node2) = C(node1,node2) - capacitance;
C(node2,node1) = C(node2,node1) - capacitance;
end
