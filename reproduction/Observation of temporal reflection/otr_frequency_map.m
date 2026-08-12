function [fout,betaConserved,info] = otr_frequency_map(fin,fromState,toState,model,p,varargin)
%OTR_FREQUENCY_MAP Momentum-conserving OFF<->ON frequency translation.
%   [FOUT,BETA,INFO] = OTR_FREQUENCY_MAP(FIN,FROM,TO,MODEL,P) finds, in
%   the first pass band of TO, the frequency satisfying
%                       beta_TO(FOUT) = beta_FROM(FIN).
%   Frequencies are in Hz. FROM/TO are 'off' or 'on', MODEL is 'ideal' or
%   'realistic'. Complex lossy beta is mapped by conserving its real part,
%   consistent with the experimental retrieval in SI Section S4.
%
%   Optional name/value settings:
%     'SearchRange' [fmin fmax] (default [1 kHz, 250 MHz])
%     'GridSize'    positive integer (default 6001)
%
%   该函数不采用固定阻抗比，而直接在 S2/S3 色散曲线的第一通带内求根，
%   因而能复现文章 Fig. 2(e,f) 的守恒动量宽带频率映射。

if nargin < 2 || isempty(fromState), fromState = 'off'; end
if nargin < 3 || isempty(toState), toState = 'on'; end
if nargin < 4 || isempty(model), model = 'ideal'; end
if nargin < 5 || isempty(p), p = otr_parameters(); end
if ~(isnumeric(fin) && isreal(fin) && ~isempty(fin) && all(isfinite(fin(:))) && all(fin(:)>0))
    error('otr:frequencyMap:Input','fin must contain positive finite frequencies in Hz.');
end
fromState = localState(fromState);
toState = localState(toState);
range = [1e3,2.5*p.fReference];
nGrid = 6001;
if mod(numel(varargin),2) ~= 0
    error('otr:frequencyMap:Options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(char(varargin{k})); value = varargin{k+1};
    switch name
        case 'searchrange', range = value;
        case 'gridsize', nGrid = value;
        otherwise, error('otr:frequencyMap:Options','Unknown option ''%s''.',name);
    end
end
if ~(isnumeric(range) && isreal(range) && numel(range)==2 && ...
        all(isfinite(range)) && range(1)>0 && range(2)>range(1))
    error('otr:frequencyMap:Range','SearchRange must be [positiveMin largerMax].');
end
if ~(isnumeric(nGrid) && isreal(nGrid) && isscalar(nGrid) && ...
        isfinite(nGrid) && nGrid==round(nGrid) && nGrid>=101)
    error('otr:frequencyMap:Grid','GridSize must be an integer >= 101.');
end

[bIn,~,passIn] = otr_bloch_dispersion(fin,fromState,model,p);
target = real(bIn(:));
if any(~passIn(:))
    error('otr:frequencyMap:InputBand','All input frequencies must lie in the FROM first pass band.');
end

fg = linspace(range(1),range(2),nGrid);
[bg,~,pg] = otr_bloch_dispersion(fg,toState,model,p);
first = localFirstBand(pg);
if isempty(first)
    error('otr:frequencyMap:NoBand','No TO-state first pass band found in SearchRange.');
end
fb = fg(first);
bb = real(bg(first));
[bb,order] = sort(bb);
fb = fb(order);
[bb,uniqueIdx] = unique(bb,'stable');
fb = fb(uniqueIdx);

foutV = NaN(size(target));
residual = NaN(size(target));
for k = 1:numel(target)
    if target(k) < bb(1) || target(k) > bb(end)
        error('otr:frequencyMap:NoSolution', ...
            'beta at fin(%d) has no TO-state solution in the first pass band.',k);
    end
    guess = interp1(bb,fb,target(k),'linear');
    il = find(fb <= guess,1,'last');
    il = max(1,min(il,numel(fb)-1));
    bracket = sort(fb(il:il+1));
    fun = @(x) real(otr_bloch_dispersion(x,toState,model,p))-target(k);
    if fun(bracket(1))*fun(bracket(2)) <= 0
        foutV(k) = fzero(fun,bracket);
    else
        foutV(k) = fminbnd(@(x) abs(fun(x)),bracket(1),bracket(2));
    end
    residual(k) = fun(foutV(k));
end

fout = reshape(foutV,size(fin));
betaConserved = bIn;
info = struct('inputFrequency',fin,'outputFrequency',fout, ...
    'ratio',fout./fin,'fromState',fromState,'toState',toState, ...
    'model',lower(char(model)),'targetBeta',bIn, ...
    'residualBeta',reshape(residual,size(fin)), ...
    'targetFirstPassBand',[fb(1),fb(end)],'searchRange',range);
end

function s = localState(s)
if ~(ischar(s) || (isstring(s) && isscalar(s)))
    error('otr:frequencyMap:State','State must be ''off'' or ''on''.');
end
s = lower(strtrim(char(s)));
if ~any(strcmp(s,{'off','on'}))
    error('otr:frequencyMap:State','State must be ''off'' or ''on''.');
end
end

function idx = localFirstBand(pass)
pass = logical(pass(:).');
start = find(pass,1,'first');
if isempty(start), idx = []; return; end
stop = find(~pass(start:end),1,'first');
if isempty(stop), stop = numel(pass)+1; else, stop = start+stop-1; end
idx = start:(stop-1);
end
