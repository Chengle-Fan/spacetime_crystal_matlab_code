function demo11_coherent_time_interface()
%DEMO11_COHERENT_TIME_INTERFACE Two-input interference at a time boundary.
%
% Counter-propagating coherent waves collide with one global time switch.
% Their relative phase controls the two outgoing temporal branches. The
% script also compares D/B- and E/B-continuous idealized jump laws to show
% why the microscopic switching mechanism must be stated explicitly.

rootDir = startup_stm();
epsBefore = 1.5^2;
epsAfter = 2.5^2;
muBefore = 1;
muAfter = 1;

[Mdb,tauDB,rhoDB] = temporal_interface_matrix_jump( ...
    epsBefore,muBefore,epsAfter,muAfter,'DB');
[Meb,tauEB,rhoEB] = temporal_interface_matrix_jump( ...
    epsBefore,muBefore,epsAfter,muAfter,'EB');

phase = linspace(-pi,pi,721);
ratioDB = abs(tauDB/rhoDB);
ratioEB = abs(tauEB/rhoEB);
inputDB = [ones(size(phase));ratioDB*exp(1i*phase)];
inputEB = [ones(size(phase));ratioEB*exp(1i*phase)];
outputDB = Mdb*inputDB;
outputEB = Meb*inputEB;

[minimumDB,idDB] = min(abs(outputDB(1,:)));
[minimumEB,idEB] = min(abs(outputEB(1,:)));
fprintf(['Coherent cancellation: DB minimum %.3e at %.3f*pi; ', ...
    'EB minimum %.3e at %.3f*pi.\n'], ...
    minimumDB,phase(idDB)/pi,minimumEB,phase(idEB)/pi);

fig = figure('Color','w','Position',[70 70 1180 780]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact');

ax1 = nexttile(tl);
bar(ax1,[real(tauDB) real(rhoDB);real(tauEB) real(rhoEB)]);
set(ax1,'XTickLabel',{'D/B continuous','E/B continuous'});
ylabel(ax1,'Interface coefficient');
title(ax1,'Same material jump, different jump law');
legend(ax1,'\tau','\rho','Location','best');
grid(ax1,'on'); box(ax1,'on');

ax2 = nexttile(tl);
plot(ax2,phase/pi,abs(outputDB(1,:)),'b-', ...
    'LineWidth',1.4);
hold(ax2,'on');
plot(ax2,phase/pi,abs(outputDB(2,:)),'r--', ...
    'LineWidth',1.4);
plot(ax2,phase(idDB)/pi,minimumDB,'ko','MarkerFaceColor','k');
xlabel(ax2,'Input relative phase /\pi');
ylabel(ax2,'Output amplitude');
title(ax2,sprintf('D/B law, |a^-/a^+|=%.3f',ratioDB));
legend(ax2,'|E^+_{out}|','|E^-_{out}|','Location','best');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
plot(ax3,phase/pi,abs(outputEB(1,:)),'b-', ...
    'LineWidth',1.4);
hold(ax3,'on');
plot(ax3,phase/pi,abs(outputEB(2,:)),'r--', ...
    'LineWidth',1.4);
plot(ax3,phase(idEB)/pi,minimumEB,'ko','MarkerFaceColor','k');
xlabel(ax3,'Input relative phase /\pi');
ylabel(ax3,'Output amplitude');
title(ax3,sprintf('E/B law, |a^-/a^+|=%.3f',ratioEB));
legend(ax3,'|E^+_{out}|','|E^-_{out}|','Location','best');
grid(ax3,'on'); box(ax3,'on');

ax4 = nexttile(tl);
plot(ax4,phase/pi, ...
    sum(abs(outputDB).^2,1)/sum(abs(inputDB).^2,1), ...
    'k-','LineWidth',1.3);
hold(ax4,'on');
plot(ax4,phase/pi, ...
    sum(abs(outputEB).^2,1)/sum(abs(inputEB).^2,1), ...
    'm--','LineWidth',1.3);
xlabel(ax4,'Input relative phase /\pi');
ylabel(ax4,'Electric-amplitude norm ratio');
title(ax4,'Modulation-mediated amplification/deamplification');
legend(ax4,'D/B law','E/B law','Location','best');
grid(ax4,'on'); box(ax4,'on');

title(tl,'Coherent wave control at a photonic time interface');
outputFile = fullfile(rootDir,'output', ...
    'demo11_coherent_time_interface.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
