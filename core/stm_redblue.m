function cmap = stm_redblue(n)
%STM_REDBLUE Blue-white-red diverging colormap.

if nargin < 1
    n = 256;
end
n1 = floor(n/2);
n2 = n - n1;
blueToWhite = [linspace(0,1,n1).', linspace(0,1,n1).', ones(n1,1)];
whiteToRed = [ones(n2,1), linspace(1,0,n2).', linspace(1,0,n2).'];
cmap = [blueToWhite; whiteToRed];
end
