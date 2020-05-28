%% FIGURE 6B: Cross-temporal correlation between execution and observation activity
% code to reproduce Figure 6B
% S.J.  Jerjian, modified April 2020

corrtype  = 'Pearson';
load('data_05132020.mat');

neurons2excl=[17 18];

% pool F5
EpochAct.popn(EpochAct.popn==4) = 3;
EpochAct.popncols(4,:) = [];
EpochAct.popnlabels = {'M1-PTN','M1-UID','F5'};
EpochAct.popn(neurons2excl) = [];
ismirror(neurons2excl) = [];

% find time of average baseline period
for c=1:length(Data)
    lcdtime(c)=mean(Data(c).events(:,strcmpi('lcdon',Data(c).eventnames)));
    cuetime(c)=mean(Data(c).events(:,strcmpi('ocue',Data(c).eventnames)));
end
[~,pos(1)] = min(abs(Data(1).times-mean(lcdtime-0.25)));
[~,pos(2)] = min(abs(Data(1).times-mean(lcdtime)));
% [~,pos(1)] = min(abs(Data(1).times-mean(lcdtime)));
% [~,pos(2)] = min(abs(Data(1).times-mean(cuetime)));

% subtract average baseline and soft-normalize
for c=1:length(Data)
    bslndata(c,:) = nanmean(Data(c).frs(pos(1):pos(2),:),1);
end
bsln = mean(bslndata,1);
for c=1:length(Data)
    Data(c).frs = bsxfun(@minus,Data(c).frs,bsln);
    Data(c).frs(:,neurons2excl) = [];
end
fr    = vertcat(Data.frs);
maxfr = max(abs(fr),[],1)+5;

figure('units','centimeters','position',[10 10 14 6],'color','w')

cmap  = flipud(cbrewer('div','RdYlGn',64));
colormap(cmap)

% one object at a time
iob = 1;
% inds = find(CONDS.obj==iob);
% frExe = Data(inds(1)).frs./maxfr;
% frObs = Data(inds(2)).frs./maxfr;
% popn = EpochAct.popn;

% both objects together
frExe_temp = Data(CONDS.hm=='m' & CONDS.gng==1); 
frObs_temp = Data(CONDS.hm=='h' & CONDS.gng==1);

frExe = [frExe_temp(1).frs frExe_temp(2).frs];
frObs = [frObs_temp(1).frs frObs_temp(2).frs]; 
% frObs = [frObs_temp(2).frs frObs_temp(1).frs]; % flip Obs

ismirror = repmat(ismirror,2,1);
popn     = repmat(EpochAct.popn,2,1);

for celltype = 1:length(EpochAct.popnlabels)
    
    ax=subplot(1,length(EpochAct.popnlabels),celltype); hold on
    
    fr1 = frExe(:,popn==celltype & ismirror);
    fr2 = frObs(:,popn==celltype & ismirror);
    
    R = corr(fr1',fr2','type',corrtype);
    R = FlipMatrixDiagonal(R);
    
    hp = imagesc(ax,Data(1).times,Data(1).times,R);

%     set(hp,'EdgeAlpha',0)
    set(ax,'clim',[-0.7 1])
    title(ax,sprintf('%s (n = %d)',EpochAct.popnlabels{celltype},size(fr1,2)/2),'interpreter','latex');
    axis(ax,[Data(1).times([1 end]) Data(1).times([1 end])])
    axis square;
    
    % plot event markers, horizontal and vertical
    for e = 1:length(Data(1).eventnames)
        cc = strcmpi(Data(1).eventnames{e},evlabs);
        plot(ones(1,2)*mean(Data(1).events(:,e)),[Data(1).times(1) Data(1).times(1)+0.2],'color',evcols(cc,:),'linewidth',2,'linestyle','-')
        plot([Data(1).times(1) Data(1).times(1)+0.2],ones(1,2)*mean(Data(1).events(:,e)),'color',evcols(cc,:),'linewidth',2,'linestyle','-')
    end
    
    if celltype==2,     xlabel('Time relative to Go [s]','interpreter','latex');
    elseif celltype==1, ylabel('Time relative to Go [s]','interpreter','latex');
    end
    
    if celltype==length(EpochAct.popnlabels)
        hc = colorbar;
        set(hc,'position',[0.92 0.24 0.02 0.56],'box','off','fontsize',9,'tickdir','out')
        hc.Ticks = [-0.7 -0.5 0 0.5 1];
        hc.TickLabels = {'','-0.5','0','0.5','1'};
        hc.TickLabelInterpreter = 'latex';
        %             title(hc,'Correlation Coefficient','fontsize',12,'interpreter','latex')
    end
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    tidyaxes(ax);
    
    Ncells(celltype,iob) = size(fr1,2)/2;
end