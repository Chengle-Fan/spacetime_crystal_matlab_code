function test_smoke()
%TEST_SMOKE Fast numerical and algebraic checks for the package.

startup_stm();

% Analytic and sampled Fourier coefficients should agree.
p = stm_fig2_parameters();
table = stpwe_sample_fourier_coefficients( ...
    @(x,t) stm_fig2_epsilon(x,t,p), -1:1, -3:3, ...
    p.Lambda,p.T,512,256);
maximumCoefficientError = 0;
for m = -1:1
    for n = -3:3
        numerical = stpwe_lookup_coefficient(table,m,n);
        analytic = stm_fig2_eps_coeff(m,n,p);
        maximumCoefficientError = max(maximumCoefficientError, ...
            abs(numerical-analytic));
    end
end
assert(maximumCoefficientError < 3e-2, ...
    'Sampled space-time Fourier coefficients failed.');

% Matching through a zero-duration intermediate medium equals one direct
% temporal interface.
[direct] = temporal_interface_matrix(1,1,9,1);
[through] = temporal_multilayer_tmm(2*pi,1,1,4,1,0,9,1);
assert(norm(direct-through,'fro') < 1e-12, ...
    'Temporal matching matrices do not compose correctly.');

% Named jump laws must preserve the quantities stated in their names.
incident = [1+0.2i;-0.3+0.1i];
[Mdb] = temporal_interface_matrix_jump(2.25,1,6.25,1,'DB');
stateBefore = temporal_directional_to_db(incident,2.25,1);
stateAfter = temporal_directional_to_db(Mdb*incident,6.25,1);
assert(norm(stateAfter-stateBefore) < 1e-12, ...
    'The D/B-continuous interface does not preserve D and B.');
[Meb,tauEB,rhoEB] = temporal_interface_matrix_jump( ...
    2.25,1,6.25,1,'EB');
afterEB = Meb*incident;
stateAfterEB = temporal_directional_to_db(afterEB,6.25,1);
assert(abs(sum(afterEB)-sum(incident)) < 1e-12 && ...
    abs(stateAfterEB(2)-stateBefore(2)) < 1e-12, ...
    'The E/B-continuous interface does not preserve E and B.');
cancelled = Meb*[1;-tauEB/rhoEB];
assert(abs(cancelled(1)) < 1e-12, ...
    'Coherent two-input cancellation failed.');

% A lossless temporal monodromy has unit determinant.
U = temporal_crystal_monodromy(1.7,[1 4],[1 1],[0.4 0.6]);
assert(abs(det(U)-1) < 1e-12, ...
    'Temporal-crystal monodromy determinant is not unity.');

% Finite-crystal helper must agree with an explicit matrix power.
kFinite = 3.1;
nPeriods = 5;
finite = temporal_finite_crystal_response(kFinite,[1 4],[1 1], ...
    [0.4 0.6],nPeriods,2,1,[1;0]);
state0 = temporal_directional_to_db([1;0],2,1);
stateExplicit = (temporal_crystal_monodromy( ...
    kFinite,[1 4],[1 1],[0.4 0.6])^nPeriods)*state0;
assert(norm(finite.finalState-stateExplicit) < 1e-11, ...
    'Finite temporal-crystal propagation failed.');

% Reversed binary PTCs support a sharply matched temporal domain-wall mode.
domainWall = temporal_domain_wall_mode( ...
    linspace(0.50,0.68,301)*2*pi,[3 1],[1 1],[0.5 0.5], ...
    [1 3],[1 1],[0.5 0.5],5,5);
assert(domainWall.bestMismatch < 3e-3, ...
    'Temporal domain-wall eigenspaces did not match.');
assert(domainWall.stateNorm(domainWall.interfaceId) >= ...
    max(domainWall.stateNorm)-1e-10, ...
    'Temporal domain-wall envelope does not peak at the interface.');

% Small-grid FHS calculation of the Rice-Mele pump must be quantized.
Nk = 21;
Np = 21;
states = complex(zeros(2,1,Nk,Np));
for ik = 1:Nk
    k = -pi+(ik-1)*2*pi/Nk;
    for ip = 1:Np
        phase = (ip-1)*2*pi/Np;
        t1 = 1+0.6*cos(phase);
        t2 = 1-0.6*cos(phase);
        mass = sin(phase);
        H = [mass,t1+t2*exp(-1i*k); ...
            t1+t2*exp(1i*k),-mass];
        [V,D] = eig(H,'vector');
        [~,order] = sort(real(D));
        states(:,1,ik,ip) = V(:,order(1));
    end
end
chern = fhs_chern_number(states);
assert(abs(chern-1) < 1e-10, 'FHS Chern-number test failed.');

% Short periodic FDTD run: a uniform plane wave should remain bounded.
Nx = 128;
L = 8;
dx = L/Nx;
x = (0:Nx-1)*dx;
dt = 0.4*dx;
n = 1.5;
k = 2*pi*4/L;
xH = x+dx/2;
E0 = exp(1i*k*x);
Hhalf0 = n*exp(1i*k*(xH+(1/n)*dt/2));
cfg.x = x;
cfg.dt = dt;
cfg.nSteps = 80;
cfg.epsFun = @(xq,tq) n^2*ones(size(xq));
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'periodic';
cfg.recordEvery = 4;
cfg.probeIndices = [1 17];
fdtd = fdtd1d_db(cfg);
relativeEnergyDrift = max(abs(fdtd.energy/fdtd.energy(1)-1));
assert(relativeEnergyDrift < 0.08, 'Uniform-medium FDTD is unstable.');
assert(fdtd.sampledCourant < 1 && size(fdtd.probeE,2) == 2, ...
    'FDTD CFL audit or probe recording failed.');

fprintf('All smoke tests passed.\n');
fprintf('Maximum sampled-coefficient error: %.3e\n', ...
    maximumCoefficientError);
fprintf('Temporal domain-wall mismatch: %.3e\n', ...
    domainWall.bestMismatch);
fprintf('FHS Chern number: %.12f\n',chern);
fprintf('Maximum short-run FDTD energy drift: %.3e\n', ...
    relativeEnergyDrift);
end
