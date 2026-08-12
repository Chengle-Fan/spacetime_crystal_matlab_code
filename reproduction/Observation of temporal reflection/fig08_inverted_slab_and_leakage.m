function results = fig08_inverted_slab_and_leakage()
%FIG08_INVERTED_SLAB_AND_LEAKAGE Reproduce SI S8/S9 qualitative phenomena.
%
% No source data were published for Extended Data Figs. 8 and 9.  This
% script therefore generates equation-driven predictions, clearly separate
% from the official-data figures.

paths=otr_setup();
Zoff=50; Zon=25;
[~,Toff,Roff]=local_interface(Zoff/Zon,'voltage');
[~,Ton,Ron]=local_interface(Zon/Zoff,'charge');
f=linspace(10,70,2401)*1e6;
taus=[17 21 25]*1e-9;
lossNpPerM=0.5*log(10)/20;       % Methods ADS A=0.5 dB/m at 100 MHz
vOn=299792458/sqrt(8.36)*(Zon/Zoff);
R=zeros(numel(taus),numel(f));
for q=1:numel(taus)
    phi=2*pi*f*(Zon/Zoff)*taus(q);
    extraDistance=vOn*taus(q);
    betaOn=2*pi*f/vOn;
    dielectricNpPerM=betaOn*0.0019/2;
    totalNpPerM=lossNpPerM.*sqrt(f/100e6)+dielectricNpPerM;
    secondLoss=exp(-2*totalNpPerM*extraDistance);
    R(q,:)=abs(Ron*Toff.*exp(-1i*phi)+ ...
        Ton*Roff.*secondLoss.*exp(1i*phi));
end

dt=.2e-9; t=(-1000:6000)*dt;
Tleak=200e-9;
signal=0.25*exp(-((t-80e-9)/(18e-9)).^2).*cos(2*pi*45e6*(t-80e-9));
leak=@(tt) .35*exp(-abs(tt)/(3e-9)).*sign(tt+eps);
raw=signal+leak(t)+leak(t-Tleak)+leak(t-2*Tleak);
shift=round(Tleak/dt);
comp=raw;
comp(shift+1:end)=raw(shift+1:end)-raw(1:end-shift);

fig=figure('Color','w','Position',[80 60 1120 700]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on'); colors=lines(numel(taus));
for q=1:numel(taus)
    plot(ax,f/1e6,R(q,:),'Color',colors(q,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('OFF duration %g ns',taus(q)*1e9));
end
xlabel(ax,'ON-state input frequency (MHz)'); ylabel(ax,'|R|');
title(ax,'(a) Inverted ON-OFF-ON slab: lossy minima do not reach zero');
legend(ax,'Location','best'); grid(ax,'on');

ax=nexttile(tl);
plot(ax,t*1e6,raw,'Color',[.65 .65 .65],'LineWidth',1.0); hold(ax,'on');
plot(ax,t*1e6,comp,'b','LineWidth',1.2);
xlabel(ax,'t (microseconds)'); ylabel(ax,'Voltage');
title(ax,'(b) Video leakage cancellation V(t)-V(t-T)');
legend(ax,'uncompensated','compensated','Location','best'); grid(ax,'on');
title(tl,'Predictions for Extended Data Figs. 8 and 9 (no raw data published)');
otr_save_figure(fig,'fig08_inverted_slab_and_leakage');

results.f=f;results.taus=taus;results.R=R;results.t=t;
results.raw=raw;results.compensated=comp;
save(fullfile(paths.output,'fig08_inverted_slab_and_leakage.mat'),'results');
end

function [M,T,R]=local_interface(rootRatio,law)
switch law
 case 'charge',T=.5*(rootRatio^2+rootRatio);R=.5*(rootRatio^2-rootRatio);
 case 'voltage',T=.5*(1+rootRatio);R=.5*(1-rootRatio);
end
M=[T R;R T];
end
