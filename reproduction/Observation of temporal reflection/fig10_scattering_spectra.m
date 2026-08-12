function results = fig10_scattering_spectra()
%FIG10_SCATTERING_SPECTRA Reproduce theoretical content of main Fig. 2a,b.
%
% The public ZIP contains no processed data for panels 2a or 2b.  Curves
% here are predictions from SI S1-S3 and the two microscopic laws S16/S17.

paths=otr_setup();p=otr_parameters();
fOff=linspace(20,70,601)*1e6;
[fOn,k]=otr_frequency_map(fOff,'off','on','ideal',p, ...
    'SearchRange',[1e5 100e6]);
[~,Zoff]=otr_bloch_dispersion(fOff,'off','ideal',p);
[~,Zon]=otr_bloch_dispersion(fOn,'on','ideal',p);
[~,Ton,Ron]=otr_temporal_interface(real(Zoff),real(Zon),'charge');
[~,Toff,Roff]=otr_temporal_interface(real(Zon),real(Zoff),'voltage');

fig=figure('Color','w','Position',[70 60 1200 690]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);plot(ax,real(k),abs(Ton),'b','LineWidth',1.5);hold(ax,'on');
plot(ax,real(k),abs(Ron),'r','LineWidth',1.5);grid(ax,'on');
ylabel(ax,'Amplitude');title(ax,'(a) OFF -> ON (charge continuous)');
legend(ax,'|T_{on}|','|R_{on}|','Location','best');
ax=nexttile(tl);plot(ax,real(k),abs(Toff),'b','LineWidth',1.5);hold(ax,'on');
plot(ax,real(k),abs(Roff),'r','LineWidth',1.5);grid(ax,'on');
ylabel(ax,'Amplitude');title(ax,'(b) ON -> OFF (voltage continuous)');
legend(ax,'|T_{off}|','|R_{off}|','Location','best');
ax=nexttile(tl);plot(ax,real(k),unwrap(angle(Ton)),'b','LineWidth',1.5);hold(ax,'on');
plot(ax,real(k),unwrap(angle(Ron)),'r','LineWidth',1.5);grid(ax,'on');
xlabel(ax,'k (rad m^{-1})');ylabel(ax,'Phase (rad)');ylim(ax,[-.2 pi+.2]);
ax=nexttile(tl);plot(ax,real(k),unwrap(angle(Toff)),'b','LineWidth',1.5);hold(ax,'on');
plot(ax,real(k),unwrap(angle(Roff)),'r','LineWidth',1.5);grid(ax,'on');
xlabel(ax,'k (rad m^{-1})');ylabel(ax,'Phase (rad)');ylim(ax,[-.2 pi+.2]);
title(tl,'Temporal scattering spectra: different microscopic boundaries');
otr_save_figure(fig,'fig10_scattering_spectra');

fprintf(['Bloch scattering at %.1f MHz OFF: Ton %.4f, Ron %.4f; ' ...
    'Toff %.4f, Roff %.4f (constant-Z target: .375, -.125, 1.5, -.5).\n'], ...
 fOff(1)/1e6,Ton(1),Ron(1),Toff(1),Roff(1));
results=struct('fOff',fOff,'fOn',fOn,'k',k,'Zoff',Zoff,'Zon',Zon, ...
 'Ton',Ton,'Ron',Ron,'Toff',Toff,'Roff',Roff);
save(fullfile(paths.output,'fig10_scattering_spectra.mat'),'results');
end
