function [templateSource,pathCleanup] = use_root_templates(requiredFunctions)
%USE_ROOT_TEMPLATES  让论文复现显式调用仓库根目录的模板求解器。
% 返回每个实际解析文件的绝对路径；调用者持有 pathCleanup 期间根目录
% 保持在 MATLAB 搜索路径首位，函数结束时精确恢复调用前的完整路径顺序。

if ~iscell(requiredFunctions) || isempty(requiredFunctions) || ...
        ~all(cellfun(@(name) ischar(name) && isrow(name) && ~isempty(name), ...
        requiredFunctions))
    error('requiredFunctions 必须是非空函数名字符元胞数组。');
end

topologyDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(topologyDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));

% 始终把根目录放到搜索路径最前，避免同名外部函数抢先被解析。
addpath(rootDirectory,'-begin');
resolvedFiles = cell(size(requiredFunctions));
for functionIndex = 1:numel(requiredFunctions)
    functionName = requiredFunctions{functionIndex};
    expectedFile = fullfile(rootDirectory,strcat(functionName,'.m'));
    if ~isfile(expectedFile)
        error('缺少根目录模板文件：%s。',expectedFile);
    end
    resolvedFile = which(functionName);
    if ~strcmp(resolvedFile,expectedFile)
        error('模板路径验证失败：%s 实际解析到 %s，而不是根目录文件 %s。', ...
            functionName,resolvedFile,expectedFile);
    end
    resolvedFiles{functionIndex} = resolvedFile;
end

templateSource = struct();
templateSource.rootDirectory = rootDirectory;
templateSource.requiredFunctions = requiredFunctions;
templateSource.resolvedFiles = resolvedFiles;
templateSource.verification = ...
    '每个求解函数均通过 which() 验证为仓库根目录的唯一模板文件';

end
