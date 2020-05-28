function [p,h,stats] = permtest(y,nulldist,alphaval,tailed,addconst)
% [p,h] = permtest(y,nulldist,alphaval,tailed)
%
% Calculates p-value for Monte Carlo simulation (permutation test, no assumption of distribution
% p-value is estimated as the number of iterations where the null distribution is (greater/smaller) than
% the test statistics divided by the total number of iterations. A constant of 1 is added by default to both sides
% so that the p-value will always be >0.
%
% INPUTS:
% y:            test statistic, scalar value
% nulldist:     null (chance) distribution of values, nx1 vector (n=number of iterations)
% [alphaval]:   desired alpha value for hypothesis test (default = 0.05)
% [tailed]:     'two','upper','lower'. tailedness of test. (default = 'two')
% [addconst]:   add constant 1 before ratio calculation (default = true) 
% OUTPUTS:
% p:        estimated p-value
% h:        hypothesis test result based on given alpha value. 1 = reject, 0 = fail to reject
% stats:    .mean   mean of null distribution
%           .ci     'confidence intervals' (percentiles) of null distribution based on tailedness
%           .test   'two','upper','lower', store value
%           .numIter size of null distribution
%
% e.g.
% 
% y = 0.69;
% nulldist = 0.4 + (0.7-0.4).*rand(10000,1); % generate chance values between 0.4 and 0.7
% [p,h] = permtest(y,nulldist,0.05,'upper'); % test is y is greater than chance
%
% @Steven Jerjian  2019/01/14
%

%% run some checks

assert(isscalar(y),'Test statistic y should be a scalar value');
assert(isvector(nulldist),'Null distribution input should be a vector');

if nargin < 5,
    addconst = false;
end
if nargin < 4,
    tailed = 'two';
end

if nargin < 3,
    alphaval = 0.05;
end

N = length(nulldist);

[p,h] = deal(nan);
stats = [];
if isnan(y)
    disp('y is NaN...'); return;
end
if sum(isnan(nulldist))/N > 0.05
    disp('More than 10% of the null distribution is NaN...'); return;
end

%% estimate significance of test statistic

% NOTE: To get rough probability that y is less than null distribution, count number of null distribution
% values that are less than y - if this is very small, then y is likely to be significantly smaller than
% value expected by chance
p0 = sum(nulldist < y); % y < chance
p1 = sum(nulldist > y); % y > chance

switch tailed
    case 'two' 
        px = min(p0,p1); % take whichever sum is smaller for two-tailed test
        stats.CI   = prctile(nulldist,100*([alphaval/2 1-alphaval/2]));
    case 'lower' 
        px = p0;
        stats.CI   = prctile(nulldist,100*(1-alphaval)); 
    case 'upper'
        px = p1;
        stats.CI   = prctile(nulldist,100*alphaval);
end

stats.numdiff = px;
stats.numIter = N;
stats.mean    = mean(nulldist);
stats.test    = tailed;

% add 1 to force min(p) = 1/(N+1) instead of zero
if addconst
    px=px+1;
    N=N+1;
end

p = px/N; 

if strcmp(tailed,'two')
    p=p*2; 
end

h = p < alphaval;

