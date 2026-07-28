function demo10_temporal_domain_wall_mode()
%DEMO10_TEMPORAL_DOMAIN_WALL_MODE Localized state at a time-domain wall.
%
% Two order-reversed binary PTCs have the same bulk spectrum. Inside their
% common momentum gap, this script matches the growing Floquet eigenstate
% of AB to the decaying eigenstate of BA. The resulting envelope peaks at
% the temporal interface.

rootDir = startup_stm();
epsA = 3;
epsB = 1;
T = 1;
durations = [0.5 0.5]*T;
nLeft = 8;
nRight = 8;
kNorm = linspace(0.45,0.72,1001);
kValues = kNorm*2*pi/T;

mode = temporal_domain_wall_mode(kValues, ...
    [epsA epsB],[1 1],durations, ...
    [epsB epsA],[1 1],durations,nLeft,nRight);

fprintf(['Temporal domain-wall mode: k/(2*pi/T)=%.6f, ', ...
    'eigenspace mismatch %.3e.\n'], ...
    mode.k/(2*pi/T),mode.bestMismatch);

fig = figure('Color','w','Position',[70 70 1180 780]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1,kNorm,mode.growthMagnitude,'b-','LineWidth',1.3);
hold(ax1,'on');
plot(ax1,kNorm,mode.decayMagnitude,'r--','LineWidth',1.3);
yline(ax1,1,'k:');
xline(ax1,mode.k/(2*pi/T),'m:','LineWidth',1.2);
xlabel(ax1,'k/(2\pi/T)');
ylabel(ax1,'|\lambda|');
title(ax1,'Common momentum gap');
legend(ax1,'AB growing branch','BA decaying branch', ...
    'Location','best');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
valid = isfinite(mode.mismatch);
semilogy(ax2,kNorm(valid),mode.mismatch(valid),'k-', ...
    'LineWidth',1.3);
hold(ax2,'on');
plot(ax2,mode.k/(2*pi/T),mode.bestMismatch,'ro', ...
    'MarkerFaceColor','r');
xlabel(ax2,'k/(2\pi/T)');
ylabel(ax2,'|det(v_{grow}^{AB},v_{decay}^{BA})|');
title(ax2,'Floquet eigenspace matching');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
semilogy(ax3,mode.cellIndex,mode.stateNorm,'o-', ...
    'Color',[0.12 0.42 0.78],'LineWidth',1.4, ...
    'MarkerFaceColor',[0.12 0.42 0.78]);
hold(ax3,'on');
xline(ax3,0,'k--');
xlabel(ax3,'Temporal-cell index');
ylabel(ax3,'||[D,B]^T|| / interface value');
title(ax3,'Temporal localization envelope');
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
leftLayers = repmat([epsA epsB],1,nLeft);
rightLayers = repmat([epsB epsA],1,nRight);
epsilonLayers = [leftLayers rightLayers];
layerIndex = -numel(leftLayers):numel(rightLayers)-1;
stairs(ax4,[layerIndex layerIndex(end)+1], ...
    [epsilonLayers epsilonLayers(end)],'LineWidth',1.3);
hold(ax4,'on');
xline(ax4,0,'k--');
xlabel(ax4,'Half-period layer index');
ylabel(ax4,'\epsilon_r');
title(ax4,'AB | BA temporal domain wall');
ylim(ax4,[0.7 3.3]);
grid(ax4,'on'); box(ax4,'on');

title(tl,'Bulk-interface matching in an ideal binary PTC');
outputFile = fullfile(rootDir,'output', ...
    'demo10_temporal_domain_wall_mode.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
