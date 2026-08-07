function animate_fdtd1d(out, options)
%ANIMATE_FDTD1D  生成 FDTD 1D D/B-Yee 仿真结果的波包时空演化动画。
%
%   animate_fdtd1d(out) 以交互模式播放电场 E(x,t) 的演化动画。
%
%   animate_fdtd1d(out, options) 通过 options 结构体自定义动画参数。
%
%   ▏  输入参数
%
%     out      - fdtd1d_db 的输出结构体
%     options  - 可选结构体，支持以下字段:
%
%       saveAs        - 保存文件名，支持 .gif, .mp4, .avi（默认: '' 仅显示）
%       fields        - 显示的场: 'E', 'H', 或 'both'（默认: 'E'）
%       viewMode      - 'line'（线图）或 'waterfall'（时空瀑布图，默认: 'line'）
%       fps           - 保存视频的帧率（默认: 30）
%       speed         - 播放速度倍率（默认: 1）
%       xlim          - x 轴范围 [xmin, xmax]（默认: 自动）
%       ylim          - 场幅值范围 [ymin, ymax]（默认: 自动对称）
%       showTime      - 是否显示时间标签（默认: true）
%       showEnergy    - 是否叠加能量曲线子图（默认: false）
%       cmap          - waterfall 模式的 colormap（默认: 'parula'）
%       title         - 动画标题（默认: 自动生成）
%       figSize       - 图形窗口大小 [w, h] 像素（默认: [800, 500]）
%       skipFrames    - 跳帧播放: 每 N 帧取一帧（默认: 1）
%
%   ▏  示例
%
%     animate_fdtd1d(out);
%     animate_fdtd1d(out, struct('saveAs','wave.gif', 'fps',15));
%     animate_fdtd1d(out, struct('viewMode','waterfall', 'showEnergy',true));
%
%   ▏  依赖: fdtd1d_db

%==========================================================================
% 参数解析
%==========================================================================

if nargin < 2, options = struct(); end

saveAs       = getOpt(options, 'saveAs', '');
fields       = lower(getOpt(options, 'fields', 'e'));
viewMode     = lower(getOpt(options, 'viewMode', 'line'));
fps          = getOpt(options, 'fps', 30);
speed        = getOpt(options, 'speed', 1);
skipFrames   = max(1, round(getOpt(options, 'skipFrames', 1)));
showTime     = logical(getOpt(options, 'showTime', true));
showEnergy   = logical(getOpt(options, 'showEnergy', false));
cmapName     = getOpt(options, 'cmap', 'parula');
plotTitle    = getOpt(options, 'title', '');
figSize      = getOpt(options, 'figSize', [800, 500]);

assert(ismember(fields, {'e','h','both'}), 'fields must be ''E'', ''H'', or ''both''.');
assert(ismember(viewMode, {'line','waterfall'}), 'viewMode must be ''line'' or ''waterfall''.');

%==========================================================================
% 提取场数据
%==========================================================================

x = out.x;
t = out.t;
nFrames = length(t);
playFrames = 1:skipFrames:nFrames;

switch fields
    case 'e'
        fieldData = out.E; fieldLabel = 'E_y'; fieldTitle = 'Electric Field';
    case 'h'
        fieldData = out.H; fieldLabel = 'H_z'; fieldTitle = 'Magnetic Field';
    case 'both'
        fieldDataE = out.E; fieldDataH = out.H;
        fieldLabel = ''; fieldTitle = 'E_y and H_z';
end

%==========================================================================
% 坐标轴范围
%==========================================================================

xlimVal = getOpt(options, 'xlim', [x(1), x(end)]);

if isfield(options, 'ylim') && ~isempty(options.ylim)
    ylimVal = options.ylim;
else
    if strcmp(fields, 'both')
        maxAbsE = max(abs(real(fieldDataE(playFrames,:))), [], 'all');
        maxAbsH = max(abs(real(fieldDataH(playFrames,:))), [], 'all');
        maxAbs = max(maxAbsE, maxAbsH);
    else
        maxAbs = max(abs(real(fieldData(playFrames,:))), [], 'all');
    end
    if maxAbs == 0, maxAbs = 1; end
    ylimVal = [-1.1*maxAbs, 1.1*maxAbs];
end

%==========================================================================
% 保存设置
%==========================================================================

doSave = ~isempty(saveAs);
if doSave
    [~, ~, ext] = fileparts(saveAs);
    ext = lower(ext);
    if strcmp(ext, '.gif')
        saveFormat = 'gif';
    elseif ismember(ext, {'.mp4', '.avi'})
        saveFormat = 'video';
    else
        error('不支持的文件格式: %s，仅支持 .gif, .mp4, .avi', ext);
    end
end

%==========================================================================
% 创建图形
%==========================================================================

figVisible = 'on';
if doSave
    % 保存模式下隐藏图形以加速
    figVisible = 'off';
end

if showEnergy
    hFig = figure('Name', 'FDTD 1D Animation', 'Color', 'w', ...
                  'Position', [100 100 figSize(1), figSize(2)*1.15], ...
                  'Visible', figVisible);
    axField = subplot(2,1,1);
    axEnergy = subplot(2,1,2);
else
    hFig = figure('Name', 'FDTD 1D Animation', 'Color', 'w', ...
                  'Position', [100 100 figSize(1), figSize(2)], ...
                  'Visible', figVisible);
    axField = gca;
end

%==========================================================================
% 初始化绘图
%==========================================================================

realE = real(fieldData);
if strcmp(fields, 'both')
    realH = real(fieldDataH);
end
frameCount = 0;

if strcmp(viewMode, 'line')
    % ---- 1D 线图 ----
    hold(axField, 'on');
    if strcmp(fields, 'both')
        hLineE = plot(axField, x, realE(playFrames(1),:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'E_y');
        hLineH = plot(axField, x, realH(playFrames(1),:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'H_z');
        legend(axField, 'Location', 'northeast');
    else
        hLine = plot(axField, x, realE(playFrames(1),:), 'b-', 'LineWidth', 1.5);
    end
    plot(axField, xlimVal, [0 0], 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlabel(axField, 'x');
    ylabel(axField, fieldLabel);
    xlim(axField, xlimVal);
    ylim(axField, ylimVal);
    grid(axField, 'on');

    if showTime
        hTimeText = text(axField, ...
            xlimVal(1)+0.02*(xlimVal(2)-xlimVal(1)), ...
            ylimVal(2)-0.08*(ylimVal(2)-ylimVal(1)), ...
            '', 'FontSize', 11, 'FontWeight', 'bold', ...
            'BackgroundColor', [1 1 1 0.7]);
    end
    hold(axField, 'off');

else
    % ---- 瀑布图: 只存播放帧，逐行累积无 NaN 间隙 ----
    nPlay = length(playFrames);
    waterfallMatrix = NaN(nPlay, length(x));

    % 初始显示: 第一播放帧（双行防单行消失）
    waterfallMatrix(1, :) = realE(playFrames(1), :);

    hImg = imagesc(axField, x, t(playFrames([1 1])), ...
                   waterfallMatrix([1 1], :));
    colormap(axField, cmapName);
    cb = colorbar(axField);
    ylabel(cb, fieldLabel);
    xlabel(axField, 'x');
    ylabel(axField, 't');
    xlim(axField, xlimVal);
    ylim(axField, [t(playFrames(1)), t(playFrames(end))]);  % 固定 y 范围，防止动画过程中 y 轴跳变
    caxis(axField, ylimVal);
    set(axField, 'YDir', 'normal');

    if ~isempty(plotTitle)
        title(axField, plotTitle);
    else
        title(axField, [fieldTitle, ' — Spacetime Evolution']);
    end

    % 时间进度线
    hold(axField, 'on');
    hProgressLine = yline(axField, t(playFrames(1)), 'r--', 'LineWidth', 1.5);
    hold(axField, 'off');
end

% --- 能量子图 ---
if showEnergy
    plot(axEnergy, t(1:playFrames(1)), out.energy(1:playFrames(1)), ...
         'k-', 'LineWidth', 1.2);
    xlabel(axEnergy, 't');
    ylabel(axEnergy, 'Energy');
    grid(axEnergy, 'on');
    hold(axEnergy, 'on');
    hEnergyMarker = plot(axEnergy, t(playFrames(1)), out.energy(playFrames(1)), ...
                         'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    hold(axEnergy, 'off');
    xlim(axEnergy, [t(1), t(end)]);
end

% 全局标题
if ~isempty(plotTitle)
    sgtitle(plotTitle);
else
    sgtitle([fieldTitle, ' — Time Evolution']);
end

%==========================================================================
% 视频写入器初始化
%==========================================================================

if doSave && strcmp(saveFormat, 'video')
    if strcmp(ext, '.mp4')
        vidProfile = 'MPEG-4';
    else
        vidProfile = 'Motion JPEG AVI';
    end
    vidWriter = VideoWriter(saveAs, vidProfile);
    vidWriter.FrameRate = fps;
    open(vidWriter);
end

%==========================================================================
% 动画主循环
%==========================================================================

tStart = tic;
targetFrameInterval = skipFrames / fps / speed;

for idx = 1:length(playFrames)
    iFrame = playFrames(idx);

    % --- 播放速度控制 (交互模式) ---
    if idx > 1 && ~doSave
        elapsed = toc(tStart);
        targetTime = (idx - 1) * targetFrameInterval;
        if elapsed < targetTime
            pause(targetTime - elapsed);
        end
    end

    frameCount = frameCount + 1;

    % --- 更新场数据 ---
    if strcmp(viewMode, 'line')
        if strcmp(fields, 'both')
            set(hLineE, 'YData', realE(iFrame, :));
            set(hLineH, 'YData', realH(iFrame, :));
        else
            set(hLine, 'YData', realE(iFrame, :));
        end
        if showTime
            set(hTimeText, 'String', sprintf('t = %.3f', t(iFrame)));
        end
    else
        % 瀑布图: 逐行累积（无间隙——waterfallMatrix 仅含播放帧）
        waterfallMatrix(idx, :) = realE(iFrame, :);
        set(hImg, 'CData', waterfallMatrix(1:idx, :));
        set(hImg, 'YData', t(playFrames(1:idx)));
        set(hProgressLine, 'Value', t(iFrame));
    end

    % --- 能量子图 ---
    if showEnergy
        set(hEnergyMarker, 'XData', t(iFrame), 'YData', out.energy(iFrame));
    end

    % 强制刷新渲染
    drawnow;

    % --- 捕获并保存帧 ---
    if doSave
        if strcmp(saveFormat, 'gif')
            frameImg = getframe(hFig);
            [imIdx, cmap] = rgb2ind(frameImg.cdata, 256, 'nodither');
            if frameCount == 1
                imwrite(imIdx, cmap, saveAs, 'gif', ...
                        'LoopCount', Inf, 'DelayTime', 1/fps);
            else
                imwrite(imIdx, cmap, saveAs, 'gif', ...
                        'WriteMode', 'append', 'DelayTime', 1/fps);
            end
        else
            writeVideo(vidWriter, getframe(hFig));
        end
    end
end

%==========================================================================
% 清理
%==========================================================================

if doSave && strcmp(saveFormat, 'video')
    close(vidWriter);
    fprintf('视频已保存至: %s (%d 帧, %d fps)\n', saveAs, frameCount, fps);
elseif doSave && strcmp(saveFormat, 'gif')
    fprintf('GIF 已保存至: %s (%d 帧, %.2f s)\n', saveAs, frameCount, frameCount/fps);
end

if ~doSave
    fprintf('动画播放完成 (%d 帧).\n', frameCount);
end

end

%==========================================================================
% 辅助: 从 options 结构体提取字段，字段不存在则返回默认值
%==========================================================================

function val = getOpt(opts, fieldName, defaultVal)
if isfield(opts, fieldName) && ~isempty(opts.(fieldName))
    val = opts.(fieldName);
else
    val = defaultVal;
end
end
