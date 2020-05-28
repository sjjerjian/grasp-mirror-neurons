%% FIGURE 7: Subspace analysis of execution + observation in movement subspace
% code to reproduce Figure 7
% Movement initiation and grasp representation in premotor and primary motor cortex mirror neurons
% Jerjian, S.J., Sahani, M., & Kraskov A.
% May 2020

% S. J. Jerjian, modified May 2020

% runPCA=1;
% iob=1;

neurons2excl = [17 18]; % neurons with too few trials (after bad EMG trials removed)

if ~runPCA   
    
    load(sprintf('Fig7_data%d.mat',iob)); % just load pre-saved projections (10000 perms)
    load('data_05132020.mat','EpochAct','Subjects')
    
    % pool F5
    EpochAct.popn(EpochAct.popn==4) = 3;
    EpochAct.popncols(4,:) = [];
    EpochAct.popnlabels = {'M1-PTN','M1-UID','F5'};
    EpochAct.popn(neurons2excl) = [];

else % run PCA from scratch
        
    load('data_05132020.mat')
    
    % pool F5
    EpochAct.popn(EpochAct.popn==4) = 3;
    EpochAct.popncols(4,:) = [];
    EpochAct.popnlabels = {'M1-PTN','M1-UID','F5'};
    
    EpochAct.popn(neurons2excl) = [];
    ismirror(neurons2excl)      = [];

    params.numPerms = 10000;  % 100 for quicker runtime, use 10000 for final figure and statistics
    params.softenNorm = 0;
    
    % select MirNs
    tempdata = Data(CONDS.gng==1);
    evs = [];

    for c=1:length(tempdata)
        
        tempdata(c).population = EpochAct.popn(ismirror)';
        tempdata(c).frs(:,neurons2excl,:) = [];
        tempdata(c).frs = tempdata(c).frs(:,ismirror,1);
        
        evs = cat(3,evs,tempdata(c).events);
    end
    
    evs = squeeze(nanmean(nanmean(evs,1),3)); % mean across sessions and conditions
    
    pos1(1) = evs(strcmp(tempdata(1).eventnames,'hpr'))-0.05;
    pos1(2) = evs(strcmp(tempdata(1).eventnames,'ho'))+0.5;

    % using cue period produces much better Exe/Obs alignment
%     pos1(1) = evs(strcmp(tempdata(1).eventnames,'go'))-0.5;
%     pos1(2) = evs(strcmp(tempdata(1).eventnames,'go'));

    [~,pos1(1)] = min(abs(tempdata(1).times-pos1(1)));
    [~,pos1(2)] = min(abs(tempdata(1).times-pos1(2)));
    
    % Exe projection
    params.numPCs = 3;
    params.axConds= iob; 
    analyzeTimes{1} = tempdata(1).times(pos1(1):pos1(2)); 
    [ProjGrasp_Exe,SummGrasp_Exe] = ExeObs_PCAbootstrap(tempdata,analyzeTimes,params);

    % Obs projection
    params.numPCs  = 3;
    params.axConds = iob+2; %[3 4] %[3]
    analyzeTimes{1} = tempdata(1).times(pos1(1):pos1(2)); 
    [ProjGrasp_Obs,SummGrasp_Obs] = ExeObs_PCAbootstrap(tempdata,analyzeTimes,params);

    clearvars -except EpochAct tempdata ProjGrasp* SummGrasp* analyzeTimes params evcols evlabs iob EMGlevel
    save(sprintf('Fig7_data%d.mat',iob))

end

%% Figure plotting

% iob = 1; % 1 - PG, 2 - WHG
obs2exe = 1; % project obs onto exe axes (1), or exe to obs (0)
showevs = 1; % show event markers on trajectories

condcols = cbrewer('qual','Set1',5);
condcols = condcols(3:4,:); 

if obs2exe
    Proj  = ProjGrasp_Exe;
    Summ  = SummGrasp_Exe;
    ccond = iob+2; % Obs projected onto Exe axes
    titletext='Execution';
    projlabs = 'EEO';
else
    Proj  = ProjGrasp_Obs;
    Summ  = SummGrasp_Obs;
    ccond = iob; % Exe projected onto Obs axes
    titletext='Observation';
    projlabs = 'OOE';
end

% Go, HPR and DO are matched for all conditions, so just take first row
evs = mean(tempdata(1).events,1);
[~,go]  = min(abs(tempdata(1).times-evs(strcmp(tempdata(1).eventnames,'go'))));
[~,hpr] = min(abs(tempdata(1).times-evs(strcmp(tempdata(1).eventnames,'hpr'))));
[~,do]  = min(abs(tempdata(1).times-evs(strcmp(tempdata(1).eventnames,'do'))));

clear hh
if length(EpochAct.popnlabels)==3
    xspace = [0.3 0.1 0.12 0.27]; xw = [0.25 0.1 0.08];
    figure('units','centimeters','position',[5 5 13 9],'color','w')

elseif length(EpochAct.popnlabels)==4
    xspace = [0.22 0.1 0.1 0.22]; xw = [0.16 0.07 0.05];
    figure('units','centimeters','position',[5 5 21 8.5],'color','w')

end

% 1. Plot trajectories

for celltype=1:length(EpochAct.popnlabels)
    Exe = Proj(iob,celltype).projAlltimes;
    Obs = Proj(iob+2,celltype).projAlltimes;
    
    if size(Exe,2)==2 
        Exe(:,3)=0; Obs(:,3) = 0;
    end
    
    ax=subplot('position',[xspace(2)+(celltype-1)*xspace(1) 0.53 xw(1) 0.42]);
    hold on; axis square;
    title(EpochAct.popnlabels{celltype},'color',EpochAct.popncols(celltype,:),'fontsize',10,'fontweight','bold','interpreter','latex')
    
    % flip PC1 for M1-PTNs, just so all Exe trajectories go in same direction 
    if obs2exe==1 
        if celltype==1 || (celltype==2 && iob==2)
            Exe(:,1)=-Exe(:,1);
            Obs(:,1)=-Obs(:,1); % flip for Obs as well
        end
    end 

    % Exe
    hh(1)=plot3(Exe(:,1),Exe(:,2),Exe(:,3),'linew',2,'color',condcols(1,:));
     if (ccond==iob+2 || (ccond==iob && celltype==length(EpochAct.popnlabels)))
         arrowh(Exe(:,1),Exe(:,2),'k',500,55); 
     end
    
    % plot event markers and labels
    evs = mean(tempdata(ccond).events,1);
    if showevs && (ccond==iob+2 || (ccond==iob && celltype==length(EpochAct.popnlabels)))
        for e = 1:length(evs)
            if ~any(strcmpi(tempdata(1).eventnames{e},{'go','nogo','hpr','do'})),continue,end
            cc = strcmp(tempdata(1).eventnames{e},evlabs);
            if isnan(evs(e)),continue,end
            [~,pos] = min(abs(tempdata(1).times-evs(e)));
            
            if strcmpi(tempdata(1).eventnames{e},'hpr')
                if celltype==3, sxy=[-0.8 -0.2];
                else, sxy=[-0.3 -0.3];
                end
            elseif strcmpi(tempdata(1).eventnames{e},'do')
                sxy=[-0.3 -0.3];
            elseif strcmpi(tempdata(1).eventnames{e},'go')
                sxy=[0.2 -0.2];
            end
            scatter3(Exe(pos,1),Exe(pos,2),Exe(pos,3),60,evcols(cc,:),'o','filled')
%             text(Exe(pos,1)+sxy(1),Exe(pos,2)+sxy(2),Exe(pos,3),upper(tempdata(1).eventnames{e}),'fontsize',12,'fontweight','bold','color',evcols(cc,:),'interpreter','latex')
        end
    end


    % Obs
    hold on; 
    hh(2) = plot3(Obs(:,1),Obs(:,2),Obs(:,3),'linew',2,'color',condcols(2,:));
    if iob==1 && ((celltype==length(EpochAct.popnlabels) && ccond==iob+2) || ccond==iob)
        arrowh(Obs(:,1),Obs(:,2),'k',500,59);
    end

    evs = mean(tempdata(2).events,1);
    if showevs && ((celltype==length(EpochAct.popnlabels) && ccond==iob+2) || ccond==iob)
        for e = 1:length(evs)
            if ~any(strcmpi(tempdata(2).eventnames{e},{'go','nogo','hpr','do'})),continue,end
            cc = strcmp(tempdata(2).eventnames{e},evlabs);
            if isnan(evs(e)),continue,end
            [~,pos] = min(abs(tempdata(2).times-evs(e)));
            scatter3(Obs(pos,1),Obs(pos,2),Obs(pos,3),60,evcols(cc,:),'o','filled')
            %text(Obs(pos,1),Obs(pos,2),Obs(pos,3),upper(tempdata(2).eventnames{e}),'fontsize',12,'fontweight','bold','color',evcols(cc,:),'interpreter','latex')

        end
    end
    
    axis([-1 1 -1 1 -1 1]*1.8)
    ax.XTick = -1:1;
    ax.YTick = -1:1;
    ax.XTickLabel = ''; ax.YTickLabel = '';
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    xlabel('PC1','interpreter','latex');

    ylabel('PC2','interpreter','latex'); zlabel('PC3','interpreter','latex');
    
    if celltype==1
        if iob==1,xt=0.8; else,xt=-1; end
        text(xt,-1.0,'Exe','color',condcols(1,:),'fontweight','bold','fontsize',12,'interpreter','latex')
        text(xt,-1.4,'Obs','color',condcols(2,:),'fontweight','bold','fontsize',12,'interpreter','latex')
        
        text(-2.5, 2,'A','fontsize',14,'fontweight','bold')
        text(-2.5, -2.5,'B','fontsize',14,'fontweight','bold')
    elseif celltype==2
        text(-2.5, -2.5,'C','fontsize',14,'fontweight','bold')
    elseif celltype==3
        text(-2.5, -2.5,'D','fontsize',14,'fontweight','bold')
    elseif celltype==4
        text(-2.5, -2.5,'E','fontsize',14,'fontweight','bold')
    end
    if celltype~=1, ax.YAxis.Visible='off'; end

    tidyaxes(ax,10);

end

%# 2. Var Explained 
epoch=1;
clear h stats
for celltype=1:length(EpochAct.popnlabels)
    ax=subplot('position',[xspace(3)+(celltype-1)*xspace(1) 0.1 xw(2) 0.3]);
    hold on; box off; %axis square;

    h(celltype,1:2)=plot(cumsum(Summ.varCaptproj{epoch}(:,[iob iob+2],celltype))*100,'linewidth',2);
    set(h(celltype,1),'color',condcols(1,:),'marker','o');
    set(h(celltype,2),'color',condcols(2,:),'marker','o');
    if obs2exe
        h(celltype,3)=plot(cumsum(SummGrasp_Obs.varCaptproj{epoch}(:,ccond,celltype))*100,'linewidth',2);
    else
        h(celltype,3)=plot(cumsum(SummGrasp_Exe.varCaptproj{epoch}(:,ccond,celltype))*100,'linewidth',2);
    end
    
    set(h(celltype,3),'color','k','marker','o','linestyle','--');
    
    [pR(celltype),~,stats(celltype)] = permtest(Summ.alignment{epoch}(ccond,celltype),...
        squeeze(Summ.alignmentR{epoch}(ccond,celltype,:,1)),0.05,'upper',1);
    fprintf('%s pR=%.4f\n',EpochAct.popnlabels{celltype},pR(celltype))

    fprintf('\n');
    
    if celltype==1
        text(1.5,50,sprintf('Exe-%s',projlabs(1)),'color',condcols(1,:),'fontweight','bold','fontsize',10,'interpreter','latex')
        text(1.5,40,sprintf('Obs-%s',projlabs(2)),'color',condcols(2,:),'fontweight','bold','fontsize',10,'interpreter','latex')
        text(1.5,30,sprintf('Obs-%s',projlabs(3)),'color','k','fontweight','bold','fontsize',10,'interpreter','latex')
    end
    
    ax.XTick = 1:params.numPCs;
    ax.YTick = 0:25:100;
    ax.YTickLabel = {'0','','\%','','100'};
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    xll=xlabel('PC','interpreter','latex'); xll.Position(2) = xll.Position(2)+1.4; 

    axis([0.5 params.numPCs+0.5 0 100])
    tidyaxes(ax,10);
end


%# 3. Alignment Index

ymax = max(max(Summ.alignment{epoch}(ccond,:)));

for celltype=1:length(EpochAct.popnlabels)
    ax=subplot('position',[xspace(4)+(celltype-1)*xspace(1) 0.1 xw(3) 0.3]);
    hold on; box off;
    
    % plot distribution
    pc = squeeze(Summ.alignmentR{epoch}(iob,celltype,:));
    nplot=min(length(pc),10000);
    scatter(-0.04+0.08*rand(nplot,1),pc(randperm(length(pc),nplot)),7,[1 1 1]*.5,'filled')

    line([-0.06 0.06],ones(1,2)*Summ.alignment{epoch}(ccond,celltype),'LineWidth',3,'color',condcols(ceil(ccond/2),:));

     ptxt = sigstars(pR(celltype));
     text(0.12,0.15,sprintf('p=%.4f%s',pR(celltype),ptxt),'color','k','verti','middle','fontweight','bold','fontsize',10,'interpreter','latex','rotation',90);
    
    axis(ax,[-0.15 0.15 0 0.5])
    ax.XTick = [];
    ax.XAxis.Visible='off';
    ax.YTick = 0:0.1:0.5;
    ax.YTickLabel = {'0','','.2','','.4',''};
    ax.YAxis.TickLabelInterpreter = 'latex';
    tidyaxes(ax,10);
end

