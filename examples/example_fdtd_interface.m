function example_fdtd_interface()
%EXAMPLE_FDTD_INTERFACE Wave splitting at a single abrupt temporal boundary.
%
% A complex Gaussian wave packet in a uniform medium encounters an abrupt
% change in permittivity at t = tSwitch. The forward and backward scattered
% amplitudes are extracted via directional decomposition and compared with
% the analytic Morgenthaler coefficients.
%
% See also FDTD1D_DB, TEMPORAL_INTERFACE_MATRIX

rootDir = stm_init();
lambda0 = 1;
k0 = 2*pi/lambda0;
nBefore = 1.5;
nAfter = 2.5;
tSwitch = 12;
tEnd = 30;

dx = lambda0/32;
x = 0:dx:60;
dt = 0.80*dx;
nSteps = ceil(tEnd/dt);
x0 = 17;
sigma = 3.5;
profile = @(xq) exp(-((xq-x0)/sigma).^2).* ...
    exp(1i*k0*(xq-x0));

E0 = profile(x);
xH = x(1:end-1)+dx/2;
vBefore = 1/nBefore;
Hhalf0 = nBefore*profile(xH + vBefore*dt/2);

cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) temporal_eps(xq,tq,tSwitch, ...
    nBefore^2,nAfter^2);
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 160;
cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;
out = fdtd1d_db(cfg);

[~, beforeId] = min(abs(out.t-(tSwitch-2*out.dt)));
[~, afterId] = min(abs(out.t-(tSwitch+2*out.dt)));
Ebefore = out.E(beforeId,:);
Eafter = out.E(afterId,:);
Hafter = out.H(afterId,:);
Eplus = 0.5*(Eafter + Hafter/nAfter);
Eminus = 0.5*(Eafter - Hafter/nAfter);
tauNumeric = norm(Eplus)/norm(Ebefore);
rhoNumeric = norm(Eminus)/norm(Ebefore);
[~,tauExact,rhoExact] = temporal_interface_matrix( ...
    nBefore^2,1,nAfter^2,1);

fprintf(['Temporal interface: |tau| exact %.5f, FDTD %.5f; ', ...
    '|rho| exact %.5f, FDTD %.5f.\n'], ...
    abs(tauExact),tauNumeric,abs(rhoExact),rhoNumeric);

fig = figure('Color','w','Position',[80 80 1180 760]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl,[1 2]);
imagesc(ax1,out.x,out.t,abs(out.E));
set(ax1,'YDir','normal');
hold(ax1,'on');
yline(ax1,tSwitch,'w--','LineWidth',1.3);
xlabel(ax1,'x/\lambda_0');
ylabel(ax1,'Time');
title(ax1,'|E(x,t)|: temporal refraction and temporal reflection');
colormap(ax1,parula(256));
colorbar(ax1);

ax2 = nexttile(tl);
plot(ax2,out.x,real(out.E(beforeId,:)),'k-', ...
    'LineWidth',1.0);
hold(ax2,'on');
plot(ax2,out.x,real(Eplus),'b-', 'LineWidth',1.0);
plot(ax2,out.x,real(Eminus),'r-', 'LineWidth',1.0);
xlabel(ax2,'x/\lambda_0');
ylabel(ax2,'Field');
title(ax2,'Directional decomposition just after switching');
legend(ax2,'Before','Forward','Backward','Location','best');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
plot(ax3,out.t,out.energy/out.energy(1),'LineWidth',1.4);
hold(ax3,'on');
xline(ax3,tSwitch,'k--');
xlabel(ax3,'Time');
ylabel(ax3,'Instantaneous energy / initial');
title(ax3,'Energy exchange with the modulation');
grid(ax3,'on'); box(ax3,'on');

title(tl,'D/B Yee-FDTD for an abrupt temporal interface');
outputFile = fullfile(rootDir,'output', ...
    'example_fdtd_interface.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end

function epsr = temporal_eps(x,t,tSwitch,epsBefore,epsAfter)
if t < tSwitch
    epsr = epsBefore*ones(size(x));
else
    epsr = epsAfter*ones(size(x));
end
end
