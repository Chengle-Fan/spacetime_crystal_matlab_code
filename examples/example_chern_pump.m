function example_chern_pump()
%EXAMPLE_CHERN_PUMP Quantized Chern number of a Thouless pump.
%
% A Rice-Mele cycle provides a compact, independently quantized validation
% of the FHS Chern number routine. The modulation phase serves as a
% synthetic second dimension, forming a 2D parameter space (k, phase).
% The computed Chern number should equal 1.
%
% See also FHS_CHERN_NUMBER

rootDir = stm_init();
Nk = 61;
Np = 61;
kGrid = -pi + (0:Nk-1)*2*pi/Nk;
phaseGrid = (0:Np-1)*2*pi/Np;
t0 = 1;
delta0 = 0.6;
mass0 = 1.0;
states = complex(zeros(2,1,Nk,Np));
gap = zeros(Nk,Np);

for ik = 1:Nk
    k = kGrid(ik);
    for ip = 1:Np
        phase = phaseGrid(ip);
        t1 = t0 + delta0*cos(phase);
        t2 = t0 - delta0*cos(phase);
        mass = mass0*sin(phase);
        offDiagonal = t1 + t2*exp(-1i*k);
        H = [mass, offDiagonal; conj(offDiagonal), -mass];
        [V,D] = eig(H,'vector');
        [eigenvalues,order] = sort(real(D));
        states(:,1,ik,ip) = V(:,order(1));
        gap(ik,ip) = eigenvalues(2)-eigenvalues(1);
    end
end

[chern,curvature] = fhs_chern_number(states);
fprintf('Rice-Mele pump FHS Chern number = %.12f.\n',chern);

zakVsPhase = zeros(1,Np);
for ip = 1:Np
    product = 1;
    for ik = 1:Nk
        ikNext = mod(ik,Nk)+1;
        overlap = states(:,1,ik,ip)'*states(:,1,ikNext,ip);
        product = product*overlap/abs(overlap);
    end
    zakVsPhase(ip) = -angle(product);
end
zakUnwrapped = unwrap(zakVsPhase);

fig = figure('Color','w','Position',[100 100 1120 470]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact');

ax1 = nexttile(tl);
imagesc(ax1,phaseGrid/(2*pi),kGrid/(2*pi),curvature);
set(ax1,'YDir','normal');
xlabel(ax1,'Modulation phase / 2\pi');
ylabel(ax1,'k / 2\pi');
title(ax1,sprintf('FHS curvature, C=%.6f',chern));
colormap(ax1,parula(256));
colorbar(ax1);

ax2 = nexttile(tl);
plot(ax2,phaseGrid/(2*pi),zakUnwrapped/(2*pi), ...
    'LineWidth',1.5);
xlabel(ax2,'Modulation phase / 2\pi');
ylabel(ax2,'Unwrapped Zak phase / 2\pi');
title(ax2,'Quantized polarization pump');
grid(ax2,'on'); box(ax2,'on');

ax3 = nexttile(tl);
imagesc(ax3,phaseGrid/(2*pi),kGrid/(2*pi),gap);
set(ax3,'YDir','normal');
xlabel(ax3,'Modulation phase / 2\pi');
ylabel(ax3,'k / 2\pi');
title(ax3,sprintf('Band gap, minimum %.3f',min(gap(:))));
colormap(ax3,parula(256));
colorbar(ax3);

title(tl,'Topological invariant in a time-modulated two-band model');
outputFile = fullfile(rootDir,'output', ...
    'example_chern_pump.png');
try
    exportgraphics(fig,outputFile,'Resolution',220);
catch
    print(fig,outputFile,'-dpng','-r220');
end
fprintf('Saved %s\n',outputFile);
end
