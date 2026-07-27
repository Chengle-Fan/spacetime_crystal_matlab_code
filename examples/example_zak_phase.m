function example_zak_phase()
%EXAMPLE_ZAK_PHASE Biorthogonal Zak phase of ST-PWE continuum bands.
%
% Computes the Zak phase (1D Wilson loop) for both a static reference band
% and a driven Floquet band of the Park-Min spacetime crystal.
%
% Key diagnostic: the minimum Wilson link magnitude must be checked before
% interpreting a single-band invariant. Values below ~0.1 indicate band
% near-degeneracy, requiring a composite-band Wilson loop instead.
%
% The Zak phase depends on the unit-cell origin. The Park-Min cell starts
% at a layer boundary rather than an inversion center.
%
% See also STPWE_TRACK_BAND, ZAK_PHASE_BIORTHOGONAL, STPWE_BZ_SEWING_MATRIX

rootDir = stm_init();
p = stm_preset_modulated_slab();
Nk = 161;
kGrid = -p.g/2 + (0:Nk-1)*p.g/Nk;
omegaSeed = 0.344*p.g*p.c0;

pStatic = p;
pStatic.modDepth = 0;
sysStatic = stpwe_build_system( ...
    @(m,n) stm_fourier_modulated_slab(m,n,pStatic), [], ...
    12, 0, p.g, p.Omega);
optionsStatic.omegaWindow = [0.25 0.70]*p.g*p.c0;
optionsStatic.maxImag = 1e-7*p.g*p.c0;
optionsStatic.minM0Weight = 0;
bandStatic = stpwe_track_band(sysStatic,kGrid,omegaSeed,optionsStatic);
[zakStatic,infoStatic] = zak_phase_biorthogonal( ...
    bandStatic.R,bandStatic.L,sysStatic.Bomega, ...
    stpwe_bz_sewing_matrix(sysStatic));

sysDriven = stpwe_build_system( ...
    @(m,n) stm_fourier_modulated_slab(m,n,p), [], ...
    10, 1, p.g, p.Omega);
optionsDriven.omegaWindow = [0.15 0.70]*p.g*p.c0;
optionsDriven.maxImag = 0.03*p.g*p.c0;
optionsDriven.minM0Weight = 0.005;
bandDriven = stpwe_track_band(sysDriven,kGrid,omegaSeed,optionsDriven);
[zakDriven,infoDriven] = zak_phase_biorthogonal( ...
    bandDriven.R,bandDriven.L,sysDriven.Bomega, ...
    stpwe_bz_sewing_matrix(sysDriven));

fprintf('Static reference Zak phase = %.6f*pi, minimum link %.3e.\n', ...
    zakStatic/pi,infoStatic.minimumLinkMagnitude);
fprintf('Driven Floquet Zak phase   = %.6f*pi, minimum link %.3e.\n', ...
    zakDriven/pi,infoDriven.minimumLinkMagnitude);
if infoDriven.minimumLinkMagnitude < 0.1
    warning(['The driven-band Wilson link is small. Treat its single-band ', ...
        'Zak phase as a convergence diagnostic, not a robust invariant.']);
end

fig = figure('Color','w','Position',[80 80 1160 730]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1,kGrid/p.g,real(bandStatic.omega)/(p.g*p.c0), ...
    'k-','LineWidth',1.4);
hold(ax1,'on');
plot(ax1,kGrid/p.g,real(bandDriven.omega)/(p.g*p.c0), ...
    'b-','LineWidth',1.2);
xlabel(ax1,'k\Lambda/(2\pi)');
ylabel(ax1,'Re(\omega)\Lambda/(2\pi c)');
title(ax1,'Tracked band');
legend(ax1,'Static','Driven','Location','best');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
plot(ax2,kGrid/p.g,imag(bandDriven.omega)/(p.g*p.c0), ...
    'r-','LineWidth',1.2);
xlabel(ax2,'k\Lambda/(2\pi)');
ylabel(ax2,'Im(\omega)\Lambda/(2\pi c)');
title(ax2,'Driven-band stability');
yline(ax2,0,'k:');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
plot(ax3,kGrid/p.g,bandDriven.m0Weight,'LineWidth',1.2);
xlabel(ax3,'k\Lambda/(2\pi)');
ylabel(ax3,'m=0 participation');
title(ax3,'Floquet spectral weight');
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
plot(ax4,kGrid/p.g,abs(infoStatic.links),'k-', ...
    'LineWidth',1.2);
hold(ax4,'on');
plot(ax4,kGrid/p.g,abs(infoDriven.links),'b-', ...
    'LineWidth',1.2);
xlabel(ax4,'Link starting at k\Lambda/(2\pi)');
ylabel(ax4,'|Wilson link|');
title(ax4,sprintf('Zak/\\pi: static %.3f, driven %.3f', ...
    zakStatic/pi,zakDriven/pi));
legend(ax4,'Static','Driven','Location','best');
grid(ax4,'on'); box(ax4,'on');

title(tl,'Biorthogonal Wilson-loop calculation for ST-PWE');
outputFile = fullfile(rootDir,'output','example_zak_phase.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
