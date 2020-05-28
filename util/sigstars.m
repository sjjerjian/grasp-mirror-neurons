function ptxt = sigstars(pval,varargin)

if nargin==1
    pval_levels = [0.001 0.01 0.05 99];
    pval_labels = {'***','**','*',''};
else
    pval_levels = varargin{1};
    pval_labels = varargin{2}; 
end

pp=1;
while pval > pval_levels(pp)
    pp = pp+1;
end
ptxt = pval_labels{pp};
