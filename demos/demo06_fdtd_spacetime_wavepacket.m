function demo06_fdtd_spacetime_wavepacket()
%DEMO06_FDTD_SPACETIME_WAVEPACKET Propagation through the Fig. 2 medium.
%
% The modulation is present throughout the simulation. The initial packet
% is intentionally simple; it decomposes into the available Floquet-Bloch
% modes, which makes frequency conversion visible in the probe spectrum.

rootDir = startup_stm();
p = stm_fig2_parameters();
cells = 36;
pointsPerCell = 40;
dx = p.Lambda/pointsPerCell;
x = 0:dx:cells*p.Lambda;
dt = 0.55*dx/p.c0;
tEnd = 3*p.T;
nSteps = ceil(tEnd/dt);

k0 = 0.175*p.g;
x0 = 8*p.Lambda;
sigma = 2.5*p.Lambda;
nEffective = sqrt(0.75*p.eps1 + 0.25*p.epsc);
profile = @(xq) exp(-((xq-x0)/sigma).^2).* ...
    exp(1i*k0*(xq-x0));
E0 = profile(x);
xH = x(1:end-1)+dx/2;
vEffective = p.c0/nEffective;
Hhalf0 = nEffective*profile(xH + vEffective*dt/2);

cfg.x = x;
cfg.dt = dt;
cfg.nSteps = nSteps;
cfg.epsFun = @(xq,tq) stm_fig2_epsilon(xq,tq,p);
cfg.muFun = @(xq,tq) ones(size(xq));
cfg.E0 = E0;
cfg.Hhalf0 = Hhalf0;
cfg.boundary = 'sponge';
cfg.spongeCells = 120;
cfg.spongeStrength = 0.08;
cfg.recordEvery = 2;
out = fdtd1d_db(cfg);

[~,probeId] = min(abs(out.x-14*p.Lambda));
probe = out.E(:,probeId);
dtRecord = mean(diff(out.t));
[Sprobe,fProbe,tProbe] = stm_stft(probe,dtRecord, ...
    128,12,512);

fig = figure('Color','w','Position',[70 70 1250 790]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl,[1 2]);
imagesc(ax1,out.x/p.Lambda,out.t/p.T,abs(out.E));
set(ax1,'YDir','normal');
xlabel(ax1,'x/\Lambda');
ylabel(ax1,'t/T');
title(ax1,'Wave-packet evolution in the space-time crystal');
colormap(ax1,parula(256));
colorbar(ax1);
hold(ax1,'on');
xline(ax1,out.x(probeId)/p.Lambda,'w:','LineWidth',1.2);

ax2 = nexttile(tl);
imagesc(ax2,tProbe/p.T,fProbe*p.Lambda/p.c0, ...
    20*log10(abs(Sprobe)/max(abs(Sprobe(:)))+1e-8));
set(ax2,'YDir','normal');
ylim(ax2,[0 1.0]);
caxis(ax2,[-60 0]);
xlabel(ax2,'t/T');
ylabel(ax2,'f\Lambda/c');
title(ax2,'Probe spectrogram (Floquet sidebands)');
colormap(ax2,parula(256));
colorbar(ax2);

ax3 = nexttile(tl);
plot(ax3,out.t/p.T,out.energy/out.energy(1),'LineWidth',1.4);
xlabel(ax3,'t/T');
ylabel(ax3,'Instantaneous energy / initial');
title(ax3,'Parametric energy exchange');
grid(ax3,'on'); box(ax3,'on');

title(tl,'D/B Yee-FDTD for the Park-Min Fig. 2 permittivity');
outputFile = fullfile(rootDir,'output', ...
    'demo06_fdtd_spacetime_wavepacket.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
