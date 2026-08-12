function results = fig04_boundary_conditions()
%FIG04_BOUNDARY_CONDITIONS Reproduce SI Fig. S3 and Eqs. (S13)-(S17).
%
% The left column is capacitor connection (charge conserving).  The right
% column is capacitor removal (voltage conserving).  The construction makes
% the microscopic distinction visible instead of assuming a universal D/B
% jump law.

otr_setup();
t = linspace(-2,5,4001);             % normalized by the pre-switch period
w1 = 2*pi;

% ON example in SI: Cafter=4*Cbefore -> frequency is halved.
on.Cbefore = 1; on.Cafter = 4;
on.w2 = w1*sqrt(on.Cbefore/on.Cafter);
[~,on.T,on.R] = local_coeff(sqrt(on.Cbefore/on.Cafter),'charge');
on.v = cos(w1*t);
ids = t >= 0;
on.v(ids) = on.T*cos(on.w2*t(ids)) + on.R*cos(-on.w2*t(ids));
on.qBase = on.Cbefore*on.v;
on.qAdded = zeros(size(t));
on.qAdded(ids) = (on.Cafter-on.Cbefore)*on.v(ids);
on.qTotal = on.qBase+on.qAdded;

% OFF example in SI: Cafter=Cbefore/4 -> frequency doubles.
off.Cbefore = 4; off.Cafter = 1;
off.w2 = w1*sqrt(off.Cbefore/off.Cafter);
[~,off.T,off.R] = local_coeff(sqrt(off.Cbefore/off.Cafter),'voltage');
off.v = cos(w1*t);
off.v(ids) = off.T*cos(off.w2*t(ids)) + off.R*cos(-off.w2*t(ids));
off.qLine = off.Cafter*off.v;
off.qDisconnected = zeros(size(t));
off.qDisconnected(~ids) = (off.Cbefore-off.Cafter)*off.v(~ids);
off.qSensed = off.qLine;             % disconnected charge is not in TL state
off.qSensed(~ids) = off.Cbefore*off.v(~ids);

fig = figure('Color','w','Position',[80 60 1120 750]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
plot(ax,t,on.v,'LineWidth',1.5); hold(ax,'on'); xline(ax,0,'k--');
ylabel(ax,'Voltage'); title(ax,'(a) Connect C_p: voltage jumps');
legend(ax,sprintf('C_2/C_1=%.0f, f_2/f_1=%.1f', ...
    on.Cafter/on.Cbefore,on.w2/w1),'Location','best'); grid(ax,'on');

ax = nexttile(tl);
plot(ax,t,off.v,'LineWidth',1.5); hold(ax,'on'); xline(ax,0,'k--');
ylabel(ax,'Voltage'); title(ax,'(b) Remove C_p: voltage is continuous');
legend(ax,sprintf('C_2/C_1=%.2f, f_2/f_1=%.1f', ...
    off.Cafter/off.Cbefore,off.w2/w1),'Location','best'); grid(ax,'on');

ax = nexttile(tl);
plot(ax,t,on.qBase,'LineWidth',1.2); hold(ax,'on');
plot(ax,t,on.qAdded,'--','LineWidth',1.2);
plot(ax,t,on.qTotal,'k','LineWidth',1.7); xline(ax,0,'k--');
xlabel(ax,'t/T_1'); ylabel(ax,'Charge (normalized)');
title(ax,'(c) C_2v(0^+)=C_1v(0^-)');
legend(ax,'base capacitor','added capacitor','total','Location','best');
grid(ax,'on');

ax = nexttile(tl);
plot(ax,t,off.qLine,'LineWidth',1.2); hold(ax,'on');
plot(ax,t,off.qDisconnected,'--','LineWidth',1.2);
plot(ax,t,off.qSensed,'k','LineWidth',1.7); xline(ax,0,'k--');
xlabel(ax,'t/T_1'); ylabel(ax,'Charge (normalized)');
title(ax,'(d) v(0^+)=v(0^-), sensed charge jumps');
legend(ax,'remaining line','removed branch','TL total','Location','best');
grid(ax,'on');
title(tl,'Microscopic temporal boundary conditions (SI S5 / Fig. S3)');
otr_save_figure(fig,'fig04_boundary_conditions');

tol = 5e-3;
onJump = on.Cafter*on.v(find(ids,1))-on.Cbefore;
offJump = off.v(find(ids,1))-1;
assert(abs(onJump)<tol,'Charge-continuity audit failed.');
assert(abs(offJump)<tol,'Voltage-continuity audit failed.');
fprintf('Boundary audit: ON charge jump %.3e; OFF voltage jump %.3e.\n', ...
    onJump,offJump);
results.t = t; results.on = on; results.off = off;
save(fullfile(otr_setup().output,'fig04_boundary_conditions.mat'),'results');
end

function [M,T,R] = local_coeff(rootRatio,law)
switch law
    case 'charge'
        T = 0.5*(rootRatio^2+rootRatio);
        R = 0.5*(rootRatio^2-rootRatio);
    case 'voltage'
        T = 0.5*(1+rootRatio);
        R = 0.5*(1-rootRatio);
end
M = [T R;R T];
end
