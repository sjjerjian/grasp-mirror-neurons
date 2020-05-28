function [Projection, Summary] = ExeObs_PCAbootstrap(Data,analyzeTimes,params)
% [Projection, Summary] = ExeObs_PCAbootstrap(Data,analyzeTimes,params)
%
% Define axes based on particular section of data (e.g. Execution grasp), then project some data onto
% these axes.
%
% The general structure of this code is based on the jPCA procedure adopted by Churchland et al. 2012
%
% Inputs:
% Data                  nx1 structure, with one condition per entry, each entry contains
%       .frs            timepoints x neurons matrix of firing rates
%       .times          1xt vector of times
%       .population     nx1 vector of neural population labels
%
% analyzeTimes
% params
%       .axConds        which conditions to use for defining axes (default all)
%       .normalize      normalize data? (default = true)
%       .softenNorm     constant normalization softening factor (default = 5)
%       .numPCs         number of PCs to use (default = 5)
%       .numPerms       number of iterations for random subspace alignment (default = 100)
% Outputs:
% Projection
%       .proj           The projection onto the top PC vectors, for each epoch
%       .times          Those times that were used (e.g. movement epoch)
%       .projAlltimes   Projection of all the data onto the PCA (which were derived from a subset of the data)
%       .allTimes       All the times (same as .times, but for full time series)

% Summary
%       .varCaptproj            Variance captured by proj
%       .alignment              Alignment index for each condition to axConds
%       .alignmentR             Alignment index for random permutations
%       .dimOrd  

%% check inputs

numConds = length(Data);
numTimes = size(Data(1).frs,1);
if ~isfield(Data(1),'times')
    for c = 1:length(Data)
        Data(c).times = 1:numTimes;
    end
end

% which condition(s) to use for defining axes, default to using all
axConds = 1:length(Data);
if exist('params', 'var') && isfield(params,'axConds')
    axConds = params.axConds;
end

numPCs = 5;
if exist('params', 'var') && isfield(params,'numPCs')
    numPCs = params.numPCs;
end

% if time base isn't provided, just use all times
if ~exist('analyzeTimes', 'var') || isempty(analyzeTimes)
    disp('analyzing all times');
    analyzeTimes = {Data(1).times};
end

% do we normalize
normalize = true;
if exist('params', 'var') && isfield(params,'normalize')
    normalize = params.normalize;
end

% do we soft normalize 
% 5 means that 5 spikes gets mapped to 0.5, infinity to 1, and zero to zero.
softenNorm = 5;
if exist('params', 'var') && isfield(params,'softenNorm')
    softenNorm = params.softenNorm;
end

numPerms = 100; % number of permutations for Monte Carlo testing, reduce for speed during testing
if exist('params', 'var') && isfield(params,'numPerms')
    numPerms = params.numPerms;
end

% number of iterations, and neuron count for jackknife/bootstrap alignment
nIter = 1;
K     = 1; % -ve --> jackknife, +ve --> bootstrap
if exist('params', 'var') && isfield(params,'numPerms2')
    nIter = params.numPerms2(1);   
    K     = params.numPerms2(2);
end

%% mask over conditions used for defining PC axes, and normalize

% define samples using time bases provided
% - PC axes will be defined using time period in the first cell of analyzeTimes
for ii=1:length(analyzeTimes)
    analyzeIndices{ii} = ismember(Data(axConds(1)).times, analyzeTimes{ii});
    if size(analyzeIndices{ii},1) == 1
        analyzeIndices{ii} = analyzeIndices{ii}';  % orientation matters for the repmat below
    end
end

% analyzeMask  = repmat(analyzeIndices{1},length(axConds),1);  % used to mask axmat
% if not all conditions have the same length time base
analyzeMask = [];
for c=axConds
    analyzeMask = [analyzeMask; analyzeIndices{1}(1:length(Data(c).times))];
end

analyzeMask = logical(analyzeMask);

axmat  = vertcat(Data(axConds).frs); % data for getting PC axes
bigmat = vertcat(Data.frs);          % all data

if normalize
    ranges       = range(bigmat);       % For each neuron, the firing rate range across all times and conditions
    normFactors  = (ranges+softenNorm);
    axmat        = bsxfun(@times, axmat, 1./normFactors);
    bigmat       = bsxfun(@times, bigmat, 1./normFactors);
end

Summary.normFactors = normFactors;
Summary.dimOrd = 'PC_cond_popn';


%% Get PCs and projections

for cc = 1:max(Data(1).population)
    
    for iter=1:nIter
        
        % extract timepoints across chosen conditions on which to define axes using PCA
        smallmat = axmat(analyzeMask,Data(1).population==cc);
        
        % for bootstrapping with neurons, nn = neurons to keep
        if K>0
            if K<=1  % proportion
                kk = floor(K*size(smallmat,2));
            else    % fixed number
                kk = K;
            end
        end
        nn = sort(randperm(size(smallmat,2),kk));

        if K<0  % jackknife
            nn = setdiff(1:size(smallmat,2),nn);
        end
        smallmat = smallmat(:,nn); 

        [PCvectors,score0] = pca(smallmat);
        
        if numPCs > size(PCvectors,2)
            disp('numPCs > dimensions of data...');
            return
        end
        
        PCvectors = PCvectors(:,1:numPCs); % cut down to the right number of PCs
        score0    = score0(:,1:numPCs);
        
        % ####
        % now loop over each condition and project onto PCvectors
        % have to subtract the mean of the original matrix (smallmat) from now on
        % (which was done by default in MATLAB pca call above)
        [~,mu] = demean(smallmat,1);
        
        for c = 1:numConds
            for ii=1:length(analyzeTimes)
                %# NOTE: if c==params.axConds && ii==1, then datmat = smallmat, score1 = score0, alind1 = 1
                
                alindR = nan(1,numPerms);
               
                try % skip if data out of range (e.g. for NoGo)
                    
                    datmat = Data(c).frs(analyzeIndices{ii},Data(1).population==cc);
                    datmat = datmat(:,nn);
                    
                    normF  = normFactors(Data(1).population==cc);
                    datmat = bsxfun(@times, datmat, 1./normF(:,nn));
                    datmat = bsxfun(@minus,datmat,mu);              % mean correct!!!
                    
                    origVar1 = sum(sum(datmat.^2));                 % original variance in data
                    
                    [score1,expl1,alind1] = subproj(datmat,PCvectors,numPCs); % projection into the low-D space.
                    
                    % compute alignment due to chance
                    
                    for i = 1:numPerms
                        
                        % permutation testing - shuffle neuron order
                        % this preserves overall variance (eigenvalues) but should change the eigenvectors
                        % i.e. just sampling random dimensions in space, up to the number of dimensions in the original data,
                        % expected cov(datmat) is the identity matrix
                        
                        randmat = datmat(:,randperm(size(datmat,2)));             % shuffle neurons
                        [scoreR,~,alindR(i)] = subproj(randmat,PCvectors,numPCs); % projection into the low-D space.
                        
                        % 1b. just generate two random orthonormal subspaces
%                         r1 = -1+2*rand(size(datmat,2),numPCs);
%                         r2 = -1+2*rand(size(datmat,2),numPCs);
%                         [Vr1,~,~] = svd(demean(r1,1)); Vr1=Vr1(:,1:numPCs);
%                         [Vr2,~,~] = svd(demean(r2,1)); Vr2=Vr2(:,1:numPCs);                      
%                         alindR(i) = trace(Vr2' * (Vr1 * Vr1') * Vr2) / numPCs;
                        
                    end
                catch
                    disp('Data out of analyzeTimes range (NoGo?) ... skipping')
                    [score1,expl1,alind1] = deal(nan);
                end
                
                %# 2. project all times (for visualization)
                datmat = Data(c).frs(:,Data(1).population==cc);
                datmat = datmat(:,nn);
                
                normF  = normFactors(Data(1).population==cc);
                datmat = bsxfun(@times, datmat, 1./normF(:,nn));
                datmat = bsxfun(@minus,datmat,mu); % mean correct
                
                [scoref] = subproj(datmat,PCvectors,numPCs); % projection of all timepoints into the low-D space
                
                %# store all results
                Projection(c,cc).proj{ii}(:,:,iter)     = score1;
                Projection(c,cc).projAlltimes(:,:,iter) = scoref;
                
                if iter==1
                    Projection(c,cc).times{ii}      = analyzeTimes{ii};
                    Projection(c,cc).Alltimes       = Data(c).times';
                end
                
                Projection(c,cc).PCvectors(:,:,iter)    = PCvectors;
                Projection(c,cc).PCmeans(:,iter)        = mu;
                
                Summary.varCaptproj{ii}(:,c,cc,iter) = expl1(:,1)/100;
                
                Summary.alignment{ii}(c,cc,iter)     = alind1;
                Summary.alignmentR{ii}(c,cc,:,iter)  = alindR;    
                
            end
        end
    end
end