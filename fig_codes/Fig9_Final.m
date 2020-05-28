%% FIGURE 9: Subspace analysis of execution, observation, and NoGo activity
% code to reproduce Figure 9
% Movement initiation and grasp representation in premotor and primary motor cortex mirror neurons
% Jerjian, S.J., Sahani, M., & Kraskov A.
% May 2020

% S. J. Jerjian, modified May 2020

% runPCA = 1;
% iob = 1;

% exclude some neurons with insufficient trials (after EMG-based discard)
neurons2excl = [17 18 181 182 183 184 286 287];

if ~runPCA
    
    load(sprintf('Fig9_data%d.mat',iob)); % just load pre-saved projections (10000 perms) 
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
        
    params.numPerms = 10000; % 100 for quick runtime

    evs = [];
    
    for c=1:length(Data)
        
        Data(c).population = EpochAct.popn';
        Data(c).frs(:,neurons2excl,:) = [];
        
        if CONDS.gng(c)==1
            evs = cat(3,evs,Data(c).events);
        end
    end
    evs = squeeze(nanmean(nanmean(evs,1),3)); % mean across sessions and conditions
       
    pos(1) = evs(strcmp(Data(1).eventnames,'go'))+0.1;
    pos(2) = evs(strcmp(Data(1).eventnames,'go'))+0.4;
    [~,pos(1)] = min(abs(Data(1).times-pos(1)));
    [~,pos(2)] = min(abs(Data(1).times-pos(2)));
    
    params.numPCs = 3;
    params.axConds= iob+2;
    analyzeTimes{1} = Data(1).times(pos(1):pos(2));
    [Proj,Summ] = ExeObs_PCAbootstrap(Data,analyzeTimes,params);
    
    clearvars -except EpochAct Data Proj Summ analyzeTimes params evcols evlabs iob EMGlevel

    save(sprintf('Fig9_data%d.mat',iob))
    
end

%% Figure plotting

% iob=1; % 1 - PG, 2 - WHG
showevs=1; % show event markers on trajectories

condcols = cbrewer('qual','Set1',5);
condcols = condcols(3:end,:);

scaling = vertcat(Proj.projAlltimes);
scaling = max(abs(scaling(:)))*2;

evs = mean(Data(1).events,1);
[~,go]  = min(abs(Data(1).times-evs(strcmp(Data(1).eventnames,'go'))));
[~,hpr] = min(abs(Data(1).times-evs(strcmp(Data(1).eventnames,'hpr'))));
[~,do]  = min(abs(Data(1).times-evs(strcmp(Data(1).eventnames,'do'))));

% only show -100 to +400ms around Go
[~,posp(1)] = min(abs(Proj(1,1).Alltimes+0.1));
[~,posp(2)] = min(abs(Proj(1,1).Alltimes-0.4));

go  = go-posp(1)+1;
hpr = hpr-posp(1)+1;
do  = do-posp(1)+1;

if length(EpochAct.popnlabels)==3
    xspace = [0.3 0.1 0.12 0.27]; xw = [0.25 0.1 0.08];
    figure('units','centimeters','position',[5 5 13 9.5],'color','w')
elseif length(EpochAct.popnlabels)==4
    xspace = [0.22 0.1 0.12 0.22]; xw = [0.16 0.07 0.05];
    figure('units','centimeters','position',[5 5 21 8.5],'color','w')

end

for celltype=1:length(EpochAct.popnlabels)
    Exe = Proj(iob,celltype).projAlltimes(posp(1):posp(2),:);
    Obs = Proj(iob+2,celltype).projAlltimes(posp(1):posp(2),:);
    Nogo= Proj(iob+4,celltype).projAlltimes(posp(1):posp(2),:);

    % just flipping for consistency of visualization
    if iob==1
        if celltype==2
            Exe(:,2) = -Exe(:,2);
            Obs(:,2) = -Obs(:,2);
            Nogo(:,2)= -Nogo(:,2);
        end
    else
        if celltype==1
            Exe(:,2) = -Exe(:,2);
            Obs(:,2) = -Obs(:,2);
            Nogo(:,2)= -Nogo(:,2);
        elseif celltype==2
            Exe(:,1) = -Exe(:,1);
            Obs(:,1) = -Obs(:,1);
            Nogo(:,1) = -Nogo(:,1);
        end
    end        

    ax=subplot('position',[xspace(2)+(celltype-1)*xspace(1) 0.53 xw(1) 0.42]);
    hold on; axis square;
    
    title(EpochAct.popnlabels{celltype},'color',EpochAct.popncols(celltype,:),'fontsize',16,'fontweight','bold','interpreter','latex')

    % Exe
    hh(1) = plot3(Exe(:,1),Exe(:,2),Exe(:,3),'linew',2,'color',condcols(1,:));
 
    evs = mean(Data(iob).events,1);
    if showevs
        for e = 1:length(evs)
            if ~any(strcmpi(Data(iob).eventnames{e},{'go','nogo','hpr','do'})),continue,end
            cc = strcmp(Data(iob).eventnames{e},evlabs);
            if isnan(evs(e)),continue,end
            [~,pos] = min(abs(Data(iob).times-evs(e)));
            pos=pos-posp(1)+1;
            if pos>size(Exe,1),continue,end
            scatter3(Exe(pos,1),Exe(pos,2),Exe(pos,3),60,evcols(cc,:),'o','filled')
        end
    end
    axis([-1 1 -1 1 -1 1]*1.25)

    % Obs
    hold on; 

    hh(2) = plot3(Obs(:,1),Obs(:,2),Obs(:,3),'linew',2,'color',condcols(2,:));
    arrowh(Obs(:,1),Obs(:,2),condcols(2,:),400,80);
    
    evs = mean(Data(iob+2).events,1);
    if showevs
        for e = 1:length(evs)
            if ~any(strcmpi(Data(iob+2).eventnames{e},{'go','nogo','hpr','do'})),continue,end
            cc = strcmp(Data(iob+2).eventnames{e},evlabs);
            if isnan(evs(e)),continue,end
            [~,pos] = min(abs(Data(iob+2).times-evs(e)));
            pos=pos-posp(1)+1;
            if pos>size(Obs,1),continue,end
            scatter3(Obs(pos,1),Obs(pos,2),Obs(pos,3),60,evcols(cc,:),'o','filled')
        end
    end

    % NoGo
    hold on
    hh(3) = plot3(Nogo(:,1),Nogo(:,2),Nogo(:,3),'linew',2,'color',condcols(3,:));

    hpr_nogo = evs(strcmpi(Data(iob+2).eventnames,'hpr'));
    [~,pos]  = min(abs(Data(iob+2).times-hpr_nogo));
    pos = pos-posp(1)+1;
    cc = strcmp('hpr',evlabs);
    scatter3(Nogo(pos,1),Nogo(pos,2),Nogo(pos,3),60,evcols(cc,:),'o','filled')

    evs = mean(Data(iob+4).events,1);
    if showevs
        for e = 1:length(evs)
            if ~any(strcmpi(Data(iob+2).eventnames{e},{'go','nogo','hpr'})),continue,end

            cc = strcmp(Data(iob+2).eventnames{e},evlabs);
            if isnan(evs(e)),continue,end

            [~,pos] = min(abs(Data(iob+2).times-evs(e)));
            pos = pos-posp(1)+1;
            if pos>size(Nogo,1),continue,end
            scatter3(Nogo(pos,1),Nogo(pos,2),Nogo(pos,3),60,evcols(cc,:),'o','filled')
        end
    end

    ax.XTick = -1:1;
    ax.YTick = -1:1;
    ax.XTickLabel = ''; ax.YTickLabel = '';
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    xlabel('PC1','interpreter','latex');
    ylabel('PC2','interpreter','latex'); zlabel('PC3','interpreter','latex');
    
    
    if celltype==1
        xt=0.7;
        text(xt,-0.2,'Exe','color',condcols(1,:),'fontweight','bold','fontsize',12,'interpreter','latex')
        text(xt,-0.5,'Obs','color',condcols(2,:),'fontweight','bold','fontsize',12,'interpreter','latex')
        text(xt,-0.8,'NoGo','color',condcols(3,:),'fontweight','bold','fontsize',12,'interpreter','latex')
        
        text(-1.8, 1.35,'A','fontsize',14,'fontweight','bold')
        text(-1.8, -1.9,'B','fontsize',14,'fontweight','bold')
    elseif celltype==2
        text(-1.8, -1.9,'C','fontsize',14,'fontweight','bold')
    elseif celltype==3
        text(-1.8, -1.9,'D','fontsize',14,'fontweight','bold')
    elseif celltype==4
        text(-1.8, -1.9,'E','fontsize',14,'fontweight','bold')
    end
    
    if celltype~=1, ax.YAxis.Visible='off'; end
        
    tidyaxes(ax,10);

end

% 2. Var Explained

epoch=1;
ccond=iob+4; % NoGo condition

clear h stats
for celltype=1:length(EpochAct.popnlabels)
    ax=subplot('position',[xspace(3)+(celltype-1)*xspace(1) 0.1 xw(2) 0.3]);
    hold on; box off; %axis square;
    
    % varExp for Exe,Obs,NoGo
    h(celltype,:)=plot(cumsum(Summ.varCaptproj{epoch}(:,[0 2 4]+iob,celltype))*100,'linewidth',2);
    
    set(h(celltype,1),'color',condcols(1,:),'marker','o','linestyle','-');
    set(h(celltype,2),'color',condcols(2,:),'marker','o','linestyle','-');
    set(h(celltype,3),'color',condcols(3,:),'marker','o','linestyle','-');
    
    % Compare alignment to random dimensions via upper-tailed permutation test
    [pR_1(celltype),~,stats(celltype)] = permtest(Summ.alignment{epoch}(ccond,celltype),...
        squeeze(Summ.alignmentR{epoch}(ccond,celltype,:,1)),0.05,'upper',1);
    
    [pR_2(celltype),~,stats(celltype)] = permtest(Summ.alignment{epoch}(1,celltype),...
        squeeze(Summ.alignmentR{epoch}(1,celltype,:,1)),0.05,'upper',1);

    fprintf('\n');
    fprintf('%s pR_NoGo = %.4f\n',EpochAct.popnlabels{celltype},pR_1(celltype))
    fprintf('%s pR_Exe = %.4f\n\n',EpochAct.popnlabels{celltype},pR_2(celltype))

    ax.XTick = 1:params.numPCs;
    ax.YTick = 0:25:100;
    ax.YTickLabel = {'0','','\%','','100'};
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    xlabel('PC','interpreter','latex')

    axis([0.5 params.numPCs+0.5 0 100])
    tidyaxes(ax,10);
end

% 3. Alignment index
for celltype=1:length(EpochAct.popnlabels)
    ax=subplot('position',[xspace(4)+(celltype-1)*xspace(1) 0.1 xw(3) 0.3]);

    hold on; box off; 

    for ccond=[iob iob+4]
        pc = squeeze(Summ.alignmentR{epoch}(ccond,celltype,:));
        scatter(-0.04+0.08*rand(size(pc)),pc,10,condcols(ceil(ccond/2),:)*0.8,'filled')
    end
    for ccond=[iob iob+4]
        line([-0.07 0.07],ones(1,2)*Summ.alignment{epoch}(ccond,celltype),'LineWidth',2.5,'color',condcols(ceil(ccond/2),:));
    end
    
    ptxt1 = sigstars(pR_1(celltype));
    ptxt2 = sigstars(pR_2(celltype));
    text(0.1,0.15,sprintf('p = %.4f%s',pR_1(celltype),ptxt1),'color',condcols(ceil(ccond/2),:),'horizo','left','verti','middle','rotation',90,'interpreter','latex','fontsize',10,'fontweight','bold')
    text(0.18,0.15,sprintf('p = %.4f%s',pR_2(celltype),ptxt2),'color',condcols(1,:),'horizo','left','verti','middle','rotation',90,'interpreter','latex','fontsize',10,'fontweight','bold')

    axis(ax,[-0.15 0.15 0 0.5])
    ax.XTick = [];
    ax.XAxis.Visible='off';
    ax.YTick = 0:0.1:0.5;
    ax.YTickLabel = {'0','','.2','','.4',''};
    ax.YAxis.TickLabelInterpreter = 'latex';
    tidyaxes(ax,10);
end