function results = fig09_tlm_design()
%FIG09_TLM_DESIGN Reproduce SI Fig. S1 transmission-line design analysis.

paths=otr_setup(); p=otr_parameters();
f=linspace(0.2,100,1801)*1e6;
[betaOff,Zoff,passOff]=otr_bloch_dispersion(f,'off','ideal',p);
[betaOn,Zon,passOn]=otr_bloch_dispersion(f,'on','ideal',p);

% Design maps use theta and the first stop-band cutoff as independent
% variables.  For a capacitive T cell, cutoff obeys A=-1 at omega=omega_c,
% which fixes b_c=2(1+cos theta_c)/sin theta_c.  b is linear in omega.
thetaDeg=linspace(1,179,260);
cutoffRatio=linspace(.02,.80,220);
[TH,WC]=meshgrid(thetaDeg*pi/180,cutoffRatio);
thetaOp=TH*.5;                 % evaluate at omega/omega0=0.5
thetaCut=TH.*WC;
bCut=2*(1+cos(thetaCut))./sin(thetaCut);
bOp=bCut.*(.5./WC);
Aon=cos(thetaOp)-.5*bOp.*sin(thetaOp);
Bon=sin(thetaOp)+.5*bOp.*cos(thetaOp)-.5*bOp;
ZonNorm=abs(Bon./sqrt(max(1-Aon.^2,eps)));
valid=abs(Aon)<1 & isfinite(ZonNorm) & bCut>0;
ZonNorm(~valid)=NaN;
Ron=.5*ZonNorm.*(ZonNorm-1);   % S16 normalized to Zoff=1

% Equal-beta translation from omega_off/omega0=0.5 to the loaded line.
Aoff=cos(.5*TH);
ratioMap=NaN(size(TH));
for row=1:size(TH,1)
  for col=1:size(TH,2)
    theta0=TH(row,col); bc=bCut(row,col); target=Aoff(row,col);
    fun=@(r) cos(theta0*.5*r)-.5*(bc*.5*r/WC(row,col)).* ...
        sin(theta0*.5*r)-target;
    rootGrid=linspace(.02,1,100);
    val=arrayfun(fun,rootGrid); id=find(val(1:end-1).*val(2:end)<=0,1);
    if ~isempty(id)
      ratioMap(row,col)=fzero(fun,rootGrid(id:id+1));
    end
  end
end

fig=figure('Color','w','Position',[60 40 1240 810]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on');
plot(ax,real(betaOff)*p.d,f/p.fReference,'b','LineWidth',1.5);
plot(ax,real(betaOn)*p.d,f/p.fReference,'r','LineWidth',1.5);
xlabel(ax,'Bloch phase beta d');ylabel(ax,'omega/omega_0');
title(ax,'(a) Dispersion, SI Eq. S2');grid(ax,'on');
legend(ax,'switch OFF','switch ON','Location','best');

ax=nexttile(tl); hold(ax,'on');
plot(ax,f(passOff)/p.fReference,real(Zoff(passOff))/p.Z0,'b','LineWidth',1.5);
plot(ax,f(passOn)/p.fReference,real(Zon(passOn))/p.Z0,'r','LineWidth',1.5);
xlabel(ax,'omega/omega_0');ylabel(ax,'Z_B/Z_0');
title(ax,'(b) Bloch impedance, SI Eq. S3');grid(ax,'on');ylim(ax,[0 1.2]);

ax=nexttile(tl);
imagesc(ax,thetaDeg,cutoffRatio,abs(Ron));axis(ax,'xy');colorbar(ax);
xlabel(ax,'theta at omega_0 (deg)');ylabel(ax,'omega_c/omega_0');
title(ax,'(c) |R_{on}| at omega/omega_0=0.5');

ax=nexttile(tl);
imagesc(ax,thetaDeg,cutoffRatio,ratioMap);axis(ax,'xy');colorbar(ax);
xlabel(ax,'theta at omega_0 (deg)');ylabel(ax,'omega_c/omega_0');
title(ax,'(d) fixed-k frequency ratio');
title(tl,'Periodic loaded transmission-line design (Extended Data Fig. 1)');
otr_save_figure(fig,'fig09_tlm_design');

[~,id50]=min(abs(f-50e6));
fprintf('S1 design at 50 MHz: Z_on/Z0=%.5f, beta_on*d=%.5f rad.\n', ...
    real(Zon(id50))/p.Z0,real(betaOn(id50))*p.d);
results=struct('frequency',f,'betaOff',betaOff,'betaOn',betaOn, ...
 'Zoff',Zoff,'Zon',Zon,'thetaDeg',thetaDeg,'cutoffRatio',cutoffRatio, ...
 'reflectionMap',Ron,'frequencyRatioMap',ratioMap);
save(fullfile(paths.output,'fig09_tlm_design.mat'),'results');
end
