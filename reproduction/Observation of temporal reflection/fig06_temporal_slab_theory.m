function results = fig06_temporal_slab_theory()
%FIG06_TEMPORAL_SLAB_THEORY Reproduce the temporal Fabry-Perot phenomena.

paths = otr_setup();
Zoff = 50; Zon = 25;
[~,Ton,Ron] = local_interface(Zon/Zoff,'charge');
[~,Toff,Roff] = local_interface(Zoff/Zon,'voltage');
f = linspace(10,100,3001)*1e6;
taus = [15 25 35]*1e-9;

R = zeros(numel(taus),numel(f));
T = R;
for q = 1:numel(taus)
    phi = 2*pi*f*(Zon/Zoff)*taus(q);
    R(q,:) = abs(Roff*Ton.*exp(-1i*phi)+Toff*Ron.*exp(1i*phi));
    T(q,:) = abs(Toff*Ton.*exp(-1i*phi)+Roff*Ron.*exp(1i*phi));
end

fig = figure('Color','w','Position',[80 70 1120 720]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,[1 2]); hold(ax,'on');
colors = lines(numel(taus));
for q=1:numel(taus)
    plot(ax,f/1e6,R(q,:),'Color',colors(q,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('tau = %g ns',taus(q)*1e9));
end
xlabel(ax,'Outer-state frequency f_1 (MHz)'); ylabel(ax,'|R_{slab}|');
title(ax,'(a) Temporal Fabry-Perot reflection zeroes'); grid(ax,'on');
legend(ax,'Location','best');

ax = nexttile(tl); hold(ax,'on');
for q=1:numel(taus)
    plot(ax,f/1e6,T(q,:),'Color',colors(q,:),'LineWidth',1.5);
end
xlabel(ax,'f_1 (MHz)'); ylabel(ax,'|T_{slab}|');
title(ax,'(b) Net forward amplitude'); grid(ax,'on');

tauSweep = linspace(10,45,1401)*1e-9;
fFixed = 28.66e6;
phi = 2*pi*fFixed*(Zon/Zoff)*tauSweep;
Rtau = abs(Roff*Ton.*exp(-1i*phi)+Toff*Ron.*exp(1i*phi));
ax = nexttile(tl);
plot(ax,tauSweep*1e9,Rtau,'k','LineWidth',1.5);
xlabel(ax,'Temporal-slab duration (ns)'); ylabel(ax,'|R_{slab}|');
title(ax,'(c) Continuous tuning at k=1.8 rad m^{-1}'); grid(ax,'on');
title(tl,'Four causal paths and time-reflection interference');
otr_save_figure(fig,'fig06_temporal_slab_theory');

pathAmplitudes = [Ton*Toff Ron*Roff Ron*Toff Ton*Roff];
fprintf(['Slab paths [forward-forward, reflection-reflection, ', ...
    'TR at first, TR at second] = [%g %g %g %g].\n'],pathAmplitudes);
assert(all(pathAmplitudes(1:2)>0) && all(pathAmplitudes(3:4)<0), ...
    'Temporal-slab path polarity audit failed.');
results.f = f; results.taus = taus; results.R = R; results.T = T;
results.tauSweep = tauSweep; results.Rtau = Rtau;
results.coefficients = struct('Ton',Ton,'Ron',Ron,'Toff',Toff,'Roff',Roff);
save(fullfile(paths.output,'fig06_temporal_slab_theory.mat'),'results');
end

function [M,T,R] = local_interface(rootRatio,law)
switch law
    case 'charge'
        T=.5*(rootRatio^2+rootRatio); R=.5*(rootRatio^2-rootRatio);
    case 'voltage'
        T=.5*(1+rootRatio); R=.5*(1-rootRatio);
end
M=[T R;R T];
end
