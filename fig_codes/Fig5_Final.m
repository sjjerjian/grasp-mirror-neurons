%% FIGURE 5: Population firing rates
% code to reproduce Figure 5
% Movement initiation and grasp representation in premotor and primary motor cortex mirror neurons
% Jerjian, S.J., Sahani, M., & Kraskov A.
% May 2020

% S. J. Jerjian, modified May 2020

load('data_05132020.mat');

% pool F5
EpochAct.popn(EpochAct.popn==4) = 3;
EpochAct.popncols(4,:) = [];
EpochAct.popnlabels = {'M1-PTN','M1-UID','F5'};

plotorder = [1 2]; % 1 - Exe, 2 - Obs, 3 - NoGo

condcols = cbrewer('qual','Set1',5);
condcols = condcols(3:4,:);

neurons2excl = [17 18]; % neurons with <10 trials in at least one object (after bad EMG trials removed)
ismirror(neurons2excl,:,:)  = [];
EpochAct.popn(neurons2excl) = [];
FSclass.class(neurons2excl,:) = [];
FSclass.sign(neurons2excl,:) = [];

% iob=1; % set object

%% Pre-processing of Data

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
end
fr    = vertcat(Data.frs);
maxfr = max(abs(fr),[],1)+5;

% select object
tempdata = Data(CONDS.obj==iob);
tempdata = tempdata(plotorder);

% PLOTTING

% for sorting and separating subpopns, initialize
order   = [];
linepos = [];

% plot heatmap
cmap     = flipud(cbrewer('div','RdBu',128));
f1=figure('color','w','units','centimeters','position',[5 5 11 9]);

for c=1:length(tempdata)
    fr    = tempdata(c).frs ./ maxfr;
    fr(:,neurons2excl) = [];
    
    times = tempdata(c).times;
    
    ax1=subplot('position',[0.07 0.55-0.40*(c-1) 0.33 0.32]);
    
    if c==plotorder(1)
        
        [~,pos(1)] = min(abs(times-mean(tempdata(c).events(:,3),1))); % Go
        [~,pos(2)] = min(abs(times-(mean(tempdata(c).events(:,6),1)+0.7))); % Ho+0.7
        
        for fs=[1 2 4 5]
            fstemp = (FSclass.class(:,iob)==fs) & EpochAct.popn==celltype;
            frtemp = fr(:,fstemp);
            if isempty(frtemp),continue,end
            
            [~,maxpos] = max(abs(frtemp(pos(1):pos(2),:)));
            [mx,temporder]  = sort(maxpos);
            fstemp = find(fstemp);
            
            order = [order; fstemp(temporder)];
            linepos = [linepos length(temporder)];
        end
        linepos=cumsum(linepos);
    end
    
    frtemp = fr(:,order);
    hold on;
    colormap(cmap)
    hp=imagesc(ax1,times,1:size(frtemp,2),frtemp');
    set(ax1,'ydir','reverse','clim',[-1 1]);
    ht=title(ax1,Subjects{plotorder(c)},'color',condcols(plotorder(c),:),'interpreter','latex');
    ht.Position(2)=ht.Position(2)+1.6;
    for i=1:length(linepos)-1
        line(ax1,times([1 end]),ones(1,2)*(linepos(i)+0.5),'color','k','linewidth',1.5)
    end
    
    % put asterisk by units shown in Fig 4
    units = [6 22 64 189]; %[6 24 66 191]
    select_units = find(ismember(order,units));
    text((times(1)-0.3)*ones(size(select_units)),select_units+0.5,'*','verti','middle','fontsize',14,'interpreter','latex')
    
    axis(ax1,[times([1 end]) 0.5 size(frtemp,2)+0.5])
    yl = .5 + [0 size(frtemp,2)*0.1]; linestyle = '-';
    
    for e=1:length(tempdata(c).eventnames)
        cc=strcmpi(tempdata(c).eventnames{e},evlabs);
        plot(ax1,ones(1,2)*nanmean(tempdata(c).events(:,e)),yl,'color',evcols(cc,:),'linewidth',1.5,'linestyle',linestyle);
    end
    
    yll=ylabel('Unit \#','interpreter','latex'); 
    yll.Position(1) = yll.Position(1)+0.5;
    ax1.YTick = [1 size(frtemp,2)];
    if c==1
        set(ax1,'xticklabel',[]);
    else
        xlabel('Time relative to Go [s]','interpreter','latex','fontsize',14);
    end
    ax1.XTick = -2:1:2;
    ax1.XAxis.TickLabelInterpreter = 'latex';
    ax1.YAxis.TickLabelInterpreter = 'latex';

    if c==length(tempdata)
        hc = colorbar;
        set(hc,'position',[0.41 0.27 0.02 0.4],'xtick',-1:1:1,'box','off','fontsize',12,'tickdir','out','ticklabelinterpreter','latex');
        title(hc,'Net Normalized Activity','fontsize',12,'rotation',90,'position',[39 52 0],'interpreter','latex')
        
        hs=sgtitle(sprintf('%s, n=%d',EpochAct.popnlabels{celltype},size(frtemp,2)));
        set(hs,'interpreter','latex');
    end
    tidyaxes(ax1);
end

%% plot FF and FS class PSTHs

events=nanmean([tempdata.events],1);
events=nanmean(reshape(events,size(events,2)/2,[]),2);

[~,evpos(1)] = min(abs(tempdata(1).times-events(strcmp(tempdata(1).eventnames,'go'))));
[~,evpos(2)] = min(abs(tempdata(1).times-events(strcmp(tempdata(1).eventnames,'hpr'))));
[~,evpos(3)] = min(abs(tempdata(1).times-events(strcmp(tempdata(1).eventnames,'do'))));

clear h

for fs=1:2 % FF & FS
    ax2=subplot('position',[0.61 0.55-0.40*(fs-1) 0.35 0.32]);

    ax2.XAxis.TickLabelInterpreter = 'latex';
    ax2.YAxis.TickLabelInterpreter = 'latex';
    

    for c=1:length(tempdata)
        fr = tempdata(c).frs ./ maxfr;
        fr(:,neurons2excl) = [];

        times = tempdata(c).times;
         
        frtemp = fr(:,FSclass.class(:,iob)==fs & EpochAct.popn==celltype);
        
        ghd_FRs(:,c,fs) = mean(frtemp(evpos,:),2);
                
        if isempty(frtemp),continue,end
        
        hold on;
        h(c,:)=plotmeanse(frtemp,2,times,condcols(plotorder(c),:),'shade');
        
        if i==1, yL = [-.3 .8];
        else,    yL = [-.5 .8];
        end
        ylim(yL)
        
        if c==length(tempdata)
            for e=1:length(tempdata(c).eventnames)
                cc=strcmpi(tempdata(c).eventnames{e},evlabs);
                plot(ax2,ones(1,2)*events(e),yL(2)*[0.9 1],'color',evcols(cc,:));
                
            end
            if fs==2
                xlabel('Time relative to Go [s]','interpreter','latex','fontsize',14);
            else
                ax2.XTickLabel = [];
            end
            text(-1,yL(1)+0.1,sprintf('%s: n=%d',FSclass.labels{fs},size(frtemp,2)),'interpreter','latex','fontsize',12)
        end
        xlim(ax2,times([1 end]))
        tidyaxes;
        ax2.YTick = -0.6:0.2:1;
        ax2.YTickLabel = {'','-0.4','','0','','0.4','','0.8',''};
    end
    ax2.ActivePositionProperty  = 'OuterPosition';

    text(-1,0.6,'Exe','color',condcols(plotorder(1),:),'interpreter','latex','fontsize',12,'fontweight','bold')
    text(-1,0.43,'Obs','color',condcols(plotorder(2),:),'interpreter','latex','fontsize',12,'fontweight','bold')

end