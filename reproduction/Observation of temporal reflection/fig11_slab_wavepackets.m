function results = fig11_slab_wavepackets()
%FIG11_SLAB_WAVEPACKETS Reproduce SI Fig. S5 and main Fig. 3e predictions.

paths=otr_setup(); p=otr_parameters();
centres=[30 35 38 41 46]*1e6;
f=linspace(5,80,6001)*1e6;
tauEff=1/(4*21.2e6); % calibrated so f_off~38 MHz maps to first slab zero
[fOn]=otr_frequency_map(f,'off','on','ideal',p,'SearchRange',[1e5 100e6]);
a=-0.1875; % Roff*Ton = Toff*Ron for Zoff/Zon=2
R=a*exp(-1i*2*pi*fOn*tauEff)+a*exp(1i*2*pi*fOn*tauEff);
fwhm=17.5e6;
sigmaF=fwhm/(2*sqrt(2*log(2)));

t=(-80:.2:80)*1e-9;
waveforms=zeros(numel(centres),numel(t));spectra=zeros(numel(centres),numel(f));
for q=1:numel(centres)
  inputSpec=exp(-.5*((f-centres(q))/sigmaF).^2);
  spectra(q,:)=abs(inputSpec.*R);
  waveforms(q,:)=real(exp(1i*2*pi*centres(q)*t).* ...
      exp(-.5*(2*pi*sigmaF*t).^2))*max(spectra(q,:));
end

fig=figure('Color','w','Position',[50 40 1300 800]);
tl=tiledlayout(fig,numel(centres),2,'TileSpacing','compact','Padding','compact');
for q=1:numel(centres)
 ax=nexttile(tl);plot(ax,t*1e9,waveforms(q,:),'k');ylabel(ax,sprintf('%g MHz',centres(q)/1e6));
 if q==1,title(ax,'Total TR wave packet (analytic envelope)');end
 if q==numel(centres),xlabel(ax,'t (ns)');end;grid(ax,'on');
 ax=nexttile(tl);plot(ax,f/1e6,spectra(q,:),'Color',[.5 0 .65]);xline(ax,38,'k--');xlim(ax,[15 60]);
 if q==1,title(ax,'Spectrum: shared reflection null near 38 MHz');end
 if q==numel(centres),xlabel(ax,'f_1 (MHz)');end;grid(ax,'on');
end
title(tl,sprintf('Temporal slab acting on 17.5 MHz-FWHM packets (tau_eff=%.2f ns)',tauEff*1e9));
otr_save_figure(fig,'fig11_slab_wavepackets');

tau=linspace(15,35,601)*1e-9;widths=(5:10)*1e-9;
integrated=zeros(numel(widths),numel(tau));fFixed=28.66e6;
for w=1:numel(widths)
 sigma=2*sqrt(2*log(2))/(2*pi*widths(w));
 fw=linspace(max(1e6,fFixed-5*sigma),min(70e6,fFixed+5*sigma),1201);
 spec=exp(-.5*((fw-fFixed)/sigma).^2);
 fow=otr_frequency_map(fw,'off','on','ideal',p,'SearchRange',[1e5 100e6]);
 for q=1:numel(tau)
  rw=a*exp(-1i*2*pi*fow*tau(q))+a*exp(1i*2*pi*fow*tau(q));
  integrated(w,q)=sqrt(trapz(fw,abs(rw.*spec).^2)/trapz(fw,abs(spec).^2));
 end
end
fig=figure('Color','w','Position',[100 100 900 500]);hold on;
for w=1:numel(widths),plot(tau*1e9,integrated(w,:),'LineWidth',1.3, ...
 'DisplayName',sprintf('FWHM %g ns',widths(w)*1e9));end
xlabel('Slab duration (ns)');ylabel('RMS reflected amplitude');grid on;
title('Fig. 3e prediction: pulse-bandwidth-averaged temporal interference');legend('Location','best');
otr_save_figure(fig,'fig11b_slab_duration_sweep');
results=struct('f',f,'centres',centres,'tauEffective',tauEff,'R',R, ...
 'spectra',spectra,'t',t,'waveforms',waveforms,'tau',tau,'widths',widths, ...
 'integratedReflection',integrated);
save(fullfile(paths.output,'fig11_slab_wavepackets.mat'),'results');
end
