function demo12_ptc_convergence_audit()
%DEMO12_PTC_CONVERGENCE_AUDIT Convergence evidence for a binary PTC.
%
% The exact 2-by-2 monodromy is used as an independent reference for the
% discontinuous-waveform temporal PWE. The saved MAT file is a compact
% example of retaining errors, cutoffs, parameters, and runtime together
% with a paper figure.

rootDir = startup_stm();
epsA = 1;
epsB = 4;
dutyA = 0.5;
T = 1;
Omega = 2*pi/T;
kNorm = linspace(0.25,1.1,61);
kValues = kNorm*Omega;
reference = temporal_crystal_bands(kValues,[epsA epsB], ...
    [1 1],[dutyA 1-dutyA]*T);
Mvalues = [3 5 9 13 19];
errors = nan(2,numel(kValues),numel(Mvalues));
runtime = zeros(size(Mvalues));

for iM = 1:numel(Mvalues)
    timer = tic;
    Mtime = Mvalues(iM);
    epsCoeff = @(m,n) double(n == 0)*temporal_binary_eps_coeff( ...
        m,epsA,epsB,dutyA);
    muCoeff = @(m,n) double(m == 0 && n == 0);
    sys = stpwe_build_system(epsCoeff,muCoeff,0,Mtime,1,Omega);

    for ik = 1:numel(kValues)
        sol = stpwe_solve_omega(sys,kValues(ik));
        folded = stpwe_fold_frequency(sol.omega,Omega);
        validIds = find(isfinite(folded) & ...
            abs(imag(folded)) < 0.45*Omega);
        used = false(size(validIds));
        for branch = 1:2
            target = reference.omegaF(branch,ik);
            cost = abs(folded(validIds)-target);
            cost(used) = inf;
            [bestCost,localId] = min(cost);
            if isfinite(bestCost)
                errors(branch,ik,iM) = bestCost/Omega;
                used(localId) = true;
            end
        end
    end
    runtime(iM) = toc(timer);
end

medianError = squeeze(median( ...
    reshape(errors,[],numel(Mvalues)),1,'omitnan'));
maximumError = squeeze(max( ...
    reshape(errors,[],numel(Mvalues)),[],1,'omitnan'));
fprintf('Temporal-PWE convergence against exact monodromy:\n');
for iM = 1:numel(Mvalues)
    fprintf('  M=%2d: median %.3e, max %.3e, %.2f s\n', ...
        Mvalues(iM),medianError(iM),maximumError(iM),runtime(iM));
end

fig = figure('Color','w','Position',[70 70 1180 780]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl);
semilogy(ax1,Mvalues,medianError,'o-','LineWidth',1.4, ...
    'MarkerFaceColor',[0.12 0.42 0.78]);
hold(ax1,'on');
semilogy(ax1,Mvalues,maximumError,'s--','LineWidth',1.4, ...
    'MarkerFaceColor',[0.82 0.18 0.18]);
xlabel(ax1,'Temporal Fourier cutoff M');
ylabel(ax1,'|\Delta\omega|/\Omega');
title(ax1,'PWE--TMM spectral error');
legend(ax1,'Median','Maximum','Location','best');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
imagesc(ax2,kNorm,Mvalues, ...
    squeeze(max(errors,[],1,'omitnan')).');
set(ax2,'YDir','normal');
xlabel(ax2,'k/\Omega');
ylabel(ax2,'Temporal Fourier cutoff M');
title(ax2,'Worst branch error at each k');
colorbar(ax2);

ax3 = nexttile(tl);
plot(ax3,kNorm,imag(reference.omegaF.')/Omega, ...
    'k-','LineWidth',1.3);
xlabel(ax3,'k/\Omega');
ylabel(ax3,'Im(\omega_F)/\Omega');
title(ax3,'Exact growth/decay reference');
yline(ax3,0,'k:');
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
plot(ax4,Mvalues,runtime,'d-','LineWidth',1.4, ...
    'MarkerFaceColor',[0.38 0.62 0.22]);
xlabel(ax4,'Temporal Fourier cutoff M');
ylabel(ax4,'Runtime (s)');
title(ax4,'Cost recorded with accuracy');
grid(ax4,'on'); box(ax4,'on');

title(tl,'Research convergence audit: discontinuous binary PTC');
outputFile = fullfile(rootDir,'output', ...
    'demo12_ptc_convergence_audit.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end

dataFile = fullfile(rootDir,'output', ...
    'demo12_ptc_convergence_data.mat');
parameters = struct('epsA',epsA,'epsB',epsB, ...
    'dutyA',dutyA,'T',T,'kNorm',kNorm);
save(dataFile,'parameters','Mvalues','errors','medianError', ...
    'maximumError','runtime');
fprintf('Saved %s and %s\n',outputFile,dataFile);
end
