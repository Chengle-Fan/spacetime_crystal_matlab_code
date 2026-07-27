function stm_selftest()
%STM_SELFTEST Fast numerical and algebraic self-test for the STM toolbox.
%
% Run this after installation to verify that the core engine, TMM, FDTD,
% and topology routines produce correct results within expected tolerances.

stm_init();

% Analytic and sampled Fourier coefficients should agree.
p = stm_preset_modulated_slab();
table = stpwe_sample_fourier_coefficients( ...
    @(x,t) stm_permittivity_modulated_slab(x,t,p), -1:1, -3:3, ...
    p.Lambda,p.T,512,256);
maximumCoefficientError = 0;
for m = -1:1
    for n = -3:3
        numerical = stpwe_lookup_coefficient(table,m,n);
        analytic = stm_fourier_modulated_slab(m,n,p);
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

% A lossless temporal monodromy has unit determinant.
U = temporal_crystal_monodromy(1.7,[1 4],[1 1],[0.4 0.6]);
assert(abs(det(U)-1) < 1e-12, ...
    'Temporal-crystal monodromy determinant is not unity.');

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
fdtd = fdtd1d_db(cfg);
relativeEnergyDrift = max(abs(fdtd.energy/fdtd.energy(1)-1));
assert(relativeEnergyDrift < 0.08, 'Uniform-medium FDTD is unstable.');

fprintf('All self-tests passed.\n');
fprintf('Maximum sampled-coefficient error: %.3e\n', ...
    maximumCoefficientError);
fprintf('FHS Chern number: %.12f\n',chern);
fprintf('Maximum short-run FDTD energy drift: %.3e\n', ...
    relativeEnergyDrift);
end
