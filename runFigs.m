% Movement initiation and grasp representation in premotor and primary motor cortex mirror neurons
% Jerjian, S.J., Sahani, M., & Kraskov A.
% May 2020

% S.J.Jerjian, May 2020

%% 
addpath(genpath('./util'))
addpath('./fig_codes')
cd ./data

%% Figure 5
celltype = 1;  % 1 - M1-PTN, 2 - M1-UID, 3 - F5
iob = 1;       % 1 - PG, 2 - WHG
Fig5_Final

% for celltype=1:3
%     for iob=1:2
%         Fig5_Final
%         close;
%     end
% end

%% Fig 6

Fig6_Final
Fig6b_Final

%% Figure 7

% these can also be set at the top of the Fig7_Final script, be careful with
% overwriting!
runPCA = 0; % rerun PCA, otherwise use saved projections
iob = 1; 

Fig7_Final

%% Figure 9

runPCA = 0; % rerun PCA, otherwise use saved projections
iob = 1; 

Fig9_Final
