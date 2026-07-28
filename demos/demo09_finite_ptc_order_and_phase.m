function demo09_finite_ptc_order_and_phase()
%DEMO09_FINITE_PTC_ORDER_AND_PHASE Finite-duration PTC observables.
%
% The AB and BA temporal unit cells have identical infinite-crystal bands
% because their monodromies are cyclic permutations. A finite experiment,
% however, also resolves the complex phase of the time-refracted and
% time-reflected outputs. That phase changes when the temporal origin is
% shifted. It is an observable diagnostic, but is not by itself a Zak
% invariant.

rootDir = startup_stm();
epsA = 3;
epsB = 1;
muA = 1;
muB = 1;
T = 1;
durations = [0.5 0.5]*T;
nPeriods = 8;
epsBackground = 2;
muBackground = 1;

kNorm = linspace(0.35,0.82,401);
kValues = kNorm*2*pi/T;
bandsAB = temporal_crystal_bands(kValues,[epsA epsB], ...
    [muA muB],durations);
bandsBA = temporal_crystal_bands(kValues,[epsB epsA], ...
    [muB muA],durations);
responseAB = temporal_finite_crystal_response(kValues, ...
    [epsA epsB],[muA muB],durations,nPeriods, ...
    epsBackground,muBackground,[1;0]);
responseBA = temporal_finite_crystal_response(kValues, ...
    [epsB epsA],[muB muA],durations,nPeriods, ...
    epsBackground,muBackground,[1;0]);

bandDifference = max(abs(sort(bandsAB.lambda,1)- ...
    sort(bandsBA.lambda,1)),[],'all');
fprintf('AB/BA maximum Floquet-multiplier difference: %.3e.\n', ...
    bandDifference);

gainAB = sqrt(abs(responseAB.forward).^2+ ...
    abs(responseAB.backward).^2);
gainBA = sqrt(abs(responseBA.forward).^2+ ...
    abs(responseBA.backward).^2);
phaseDifference = angle(exp(1i*(responseAB.relativePhase- ...
    responseBA.relativePhase)));

fig = figure('Color','w','Position',[70 70 1180 780]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1,kNorm,real(bandsAB.halfTrace),'k-','LineWidth',1.4);
hold(ax1,'on');
yline(ax1,1,'r:');
yline(ax1,-1,'r:');
xlabel(ax1,'k/(2\pi/T)');
ylabel(ax1,'Tr(U)/2');
title(ax1,'Common infinite-crystal band discriminant');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
semilogy(ax2,kNorm,gainAB,'b-','LineWidth',1.4);
hold(ax2,'on');
semilogy(ax2,kNorm,gainBA,'r--','LineWidth',1.2);
xlabel(ax2,'k/(2\pi/T)');
ylabel(ax2,'Output electric-amplitude norm');
title(ax2,sprintf('Finite response after %d periods',nPeriods));
legend(ax2,'AB unit cell','BA unit cell','Location','best');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
plot(ax3,kNorm,responseAB.relativePhase/pi,'b-', ...
    'LineWidth',1.2);
hold(ax3,'on');
plot(ax3,kNorm,responseBA.relativePhase/pi,'r--', ...
    'LineWidth',1.2);
xlabel(ax3,'k/(2\pi/T)');
ylabel(ax3,'arg(E^-/E^+)/\pi');
title(ax3,'Time-reflected/refracted relative phase');
legend(ax3,'AB unit cell','BA unit cell','Location','best');
ylim(ax3,[-1 1]);
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
plot(ax4,kNorm,phaseDifference/pi,'m-','LineWidth',1.3);
xlabel(ax4,'k/(2\pi/T)');
ylabel(ax4,'Wrapped phase difference /\pi');
title(ax4,'Order-sensitive finite-crystal observable');
ylim(ax4,[-1 1]);
grid(ax4,'on'); box(ax4,'on');

title(tl,'Finite photonic time crystal: same bands, different phase');
outputFile = fullfile(rootDir,'output', ...
    'demo09_finite_ptc_order_and_phase.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
