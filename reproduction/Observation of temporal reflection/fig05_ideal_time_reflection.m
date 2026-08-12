function results = fig05_ideal_time_reflection()
%FIG05_IDEAL_TIME_REFLECTION Show temporal reversal and conserved momentum.

paths = otr_setup();
opts.direction = 'on';
opts.tSwitch = 45e-9;
opts.t = linspace(0,150e-9,751);
out = otr_synthetic_wavepacket(opts);

fig = figure('Color','w','Position',[70 70 1120 720]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,[1 2]);
imagesc(ax,out.x,out.t*1e9,abs(out.field)); set(ax,'YDir','normal');
hold(ax,'on'); yline(ax,out.tSwitch*1e9,'w--','LineWidth',1.2);
xlabel(ax,'z (m)'); ylabel(ax,'t (ns)');
title(ax,'(a) A spatially uniform switch creates forward and backward waves');
colorbar(ax); colormap(ax,parula(256));

[~,before] = min(abs(out.t-(out.tSwitch-2e-9)));
[~,after] = min(abs(out.t-(out.tSwitch+32e-9)));
ax = nexttile(tl);
plot(ax,out.x,real(out.field(before,:)),'k','LineWidth',1.2); hold(ax,'on');
plot(ax,out.x,real(out.forward(after,:)),'b','LineWidth',1.2);
plot(ax,out.x,real(out.backward(after,:)),'r','LineWidth',1.2);
xlabel(ax,'z (m)'); ylabel(ax,'Re(V)'); grid(ax,'on');
title(ax,'(b) Directional components');
legend(ax,'before','time-refracted','time-reflected','Location','best');

ax = nexttile(tl);
AkBefore = abs(fftshift(fft(out.initial)));
k = fftshift(out.k);
[~,atSwitch] = min(abs(out.t-out.tSwitch));
AkAfter = abs(fftshift(fft(out.field(atSwitch,:))));
keep = abs(k-out.k0) <= max(8,0.8*out.k0);
plot(ax,k(keep),AkBefore(keep)/max(AkBefore),'k','LineWidth',1.4); hold(ax,'on');
plot(ax,k(keep),AkAfter(keep)/max(AkAfter),'g--','LineWidth',1.0);
xlabel(ax,'k (rad m^{-1})'); ylabel(ax,'Normalized spectrum'); grid(ax,'on');
title(ax,sprintf('(c) k conserved; f_2/f_1 = %.3f',out.fRatio));
legend(ax,'before switch','after switch','Location','best');
title(tl,sprintf('Ideal photonic time interface: T=%.3f, R=%.3f', ...
    out.Tfirst,out.Rfirst));
otr_save_figure(fig,'fig05_ideal_time_reflection');

results = out;
save(fullfile(paths.output,'fig05_ideal_time_reflection.mat'),'results');
end
