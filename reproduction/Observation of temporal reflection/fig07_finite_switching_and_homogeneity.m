function results = fig07_finite_switching_and_homogeneity()
%FIG07_FINITE_SWITCHING_AND_HOMOGENEITY Reproduce SI Figs. S6 and S7 trends.
%
% The finite-rise curve uses a first-order Fourier-overlap estimate for a
% linear frequency ramp, |sinc(Delta omega*tau/2)|.  It is an
% explicitly labelled analytic switching-time study, not unpublished ADS
% data.  The synchronization panel evaluates phase spread across 30 cells.

paths = otr_setup();
f0 = 50e6;
ratio = 0.55;
riseNs = linspace(0,16,321);
deltaW = 2*pi*f0*(1-ratio);
x = deltaW*(riseNs*1e-9)/2;
finiteFactor = ones(size(x));
finiteFactor(x~=0) = abs(sin(x(x~=0))./x(x~=0));

% The published measurement reports ~90% at 3 ns.  Plot it as a separate
% experimental anchor rather than tuning the analytic curve to pass it.
publishedRise = 3;
publishedAmplitude = 0.90;

N = 30;
mainCellDelay = 0.208/(299792458/sqrt(8.36));
controlTotalDelay = (N-1)*mainCellDelay/80;
cellDelay = controlTotalDelay/(N-1);
cell = (0:N-1)-(N-1)/2;
frequency = linspace(10,100,901)*1e6;
coherence = zeros(size(frequency));
for n=1:numel(frequency)
    coherence(n)=abs(mean(exp(-1i*2*pi*frequency(n)*cell*cellDelay)));
end

fig=figure('Color','w','Position',[90 80 1100 450]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);
plot(ax,riseNs,finiteFactor,'LineWidth',1.6); hold(ax,'on');
plot(ax,publishedRise,publishedAmplitude,'ro','MarkerFaceColor','r', ...
    'DisplayName','paper: 3 ns, about 90%');
xlabel(ax,'10%-90% switching time (ns)'); ylabel(ax,'TR peak / abrupt limit');
title(ax,'(a) Finite switching suppresses time reflection'); grid(ax,'on');
legend(ax,'Fourier-overlap estimate','paper anchor','Location','best'); ylim(ax,[0 1.05]);

ax=nexttile(tl);
plot(ax,frequency/1e6,coherence,'k','LineWidth',1.6); hold(ax,'on');
yline(ax,1,'r--'); xlabel(ax,'Frequency (MHz)');
ylabel(ax,'Synchronized-array coherence');
title(ax,'(b) Control line 80x faster than signal line'); grid(ax,'on');
ylim(ax,[0.95 1.001]);
title(tl,'Finite rise time and spatial homogeneity audits');
otr_save_figure(fig,'fig07_finite_switching_and_homogeneity');

fprintf('Control delay across sample %.3f ns; coherence at 100 MHz %.6f.\n', ...
    controlTotalDelay*1e9,coherence(end));
results.riseNs=riseNs; results.finiteFactor=finiteFactor;
results.publishedAnchor=[publishedRise publishedAmplitude];
results.frequency=frequency; results.coherence=coherence;
results.controlTotalDelay=controlTotalDelay;
save(fullfile(paths.output,'fig07_finite_switching_and_homogeneity.mat'),'results');
end
