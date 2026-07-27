function table = stpwe_sample_fourier_coefficients(materialFun, ...
    mOrders, nOrders, Lambda, T, Nx, Nt)
%STPWE_SAMPLE_FOURIER_COEFFICIENTS Numerical coefficients for any p(x,t).
%
% materialFun(X,TIME) must return the material parameter on compatible
% arrays. Uniform endpoint-free samples are used. The coefficient convention
% is p(x,t)=sum p_mn exp(i*n*g*x-i*m*Omega*t).

if nargin < 6 || isempty(Nx)
    Nx = max(128,8*numel(nOrders));
end
if nargin < 7 || isempty(Nt)
    Nt = max(128,8*numel(mOrders));
end

x = (0:Nx-1)*Lambda/Nx;
t = (0:Nt-1)*T/Nt;
[X,Time] = meshgrid(x,t);
values = materialFun(X,Time);
if ~isequal(size(values),size(X))
    error('materialFun must return an Nt-by-Nx array.');
end

g = 2*pi/Lambda;
Omega = 2*pi/T;
coeff = complex(zeros(numel(mOrders),numel(nOrders)));
for im = 1:numel(mOrders)
    m = mOrders(im);
    for in = 1:numel(nOrders)
        n = nOrders(in);
        kernel = exp(-1i*n*g*X + 1i*m*Omega*Time);
        coeff(im,in) = mean(values.*kernel,'all');
    end
end

table.mOrders = mOrders(:).';
table.nOrders = nOrders(:).';
table.values = coeff;
table.Lambda = Lambda;
table.T = T;
table.Nx = Nx;
table.Nt = Nt;
end
