function filename = otr_save_figure(fig, basename, resolution)
%OTR_SAVE_FIGURE Export one reproduction figure with a portable fallback.

if nargin < 3 || isempty(resolution), resolution = 220; end
paths = otr_setup();
filename = fullfile(paths.output,[basename '.png']);
try
    exportgraphics(fig,filename,'Resolution',resolution);
catch
    print(fig,filename,'-dpng',sprintf('-r%d',resolution));
end
fprintf('Saved: %s\n',filename);
end
