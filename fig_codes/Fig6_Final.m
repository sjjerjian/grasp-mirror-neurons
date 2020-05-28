%% FIGURE 6: Correlation between execution and observation activity
% code to reproduce Figure 6A and Figure 6 Supplement
% Movement initiation and grasp representation in premotor and primary motor cortex mirror neurons
% Jerjian, S.J., Sahani, M., & Kraskov A.
% May 2020

% S. J. Jerjian, modified May 2020

load('data_05132020.mat','EpochAct','ismirror','CONDS'); % revision 2 (May 2020)

% ismirror: logical variable for classified MNs
% EpochAct contains all relevant information
% .EpochFRs - units * epochs * conditions
% .CONDS structure contains condition labels (object, h/m, go/nogo)

% pool F5
EpochAct.popn(EpochAct.popn==4) = 3;
EpochAct.popncols(4,:) = [];
EpochAct.popnlabels = {'M1-PTN','M1-UID','F5'};

% corrtype = 'Spearman';
corrtype  = 'Pearson';
plotsupps = 1;

% baseline-correct and soft-normalize firing rates
NormFRs = EpochAct.EpochFRs - nanmean(EpochAct.EpochFRs(:,1,:),3);
NormFRs = NormFRs ./ (max(max(abs(NormFRs),[],2),[],3)+5);

neurons2excl = [17 18]; % neurons with <10 trials in at least one object (after bad EMG trials removed)
NormFRs(neurons2excl,:,:)   = [];
ismirror(neurons2excl,:,:)  = [];
EpochAct.popn(neurons2excl) = [];

EpochAct.labs(3:7) = {'ObjCue','E React','L React','E Reach','L Reach'};
%% line plot

rng(44)

figure('units','centimeters','position',[10 10 13 9],'color','w'); clf;
ax1  = gca; hold on;

clear R P R2 P2 C Rshuff
epochs = 1:8; % epochs to show
niter  = 1000; % number of shuffle iterations

plotInset = 1;
xpos   = [0.28 0.6 0.8];

Rshuff = nan(niter,length(epochs),length(EpochAct.popnlabels));

i=0;
for e=1:length(epochs)
    for celltype=1:length(EpochAct.popnlabels)
        selc   = EpochAct.popn==celltype & ismirror;
        Nsel(celltype) = nnz(selc);
        
        exe1 = squeeze(NormFRs(selc,epochs(e),CONDS.hm=='m'&CONDS.gng==1));
        obs1 = squeeze(NormFRs(selc,epochs(e),CONDS.hm=='h'&CONDS.gng==1));   

        exe=exe1(:); obs=obs1(:); 
        
        [r,p]   = corr(exe,obs,'rows','complete','type',corrtype);
         
        if ismember(e,[3 6 8]) && celltype==1 && plotInset
           i=i+1;
           ax_inset=axes('position',[xpos(i) 0.26 0.15 0.15]);
           hold on
           scatter(ax_inset,exe1(:,1),obs1(:,1),12,'r','filled'); % PG
           scatter(ax_inset,exe1(:,2),obs1(:,2),12,'b','filled'); % WHG
           
           if e==3
               text(0.7,-0.1,'PG','color','r','interpreter','latex','fontsize',10)
               text(0.7,-0.6,'WHG','color','b','interpreter','latex','fontsize',10)
           end

           ptxt = sprintf('p = %.g',p);
           text(0,1.45,{sprintf('r = %.2f',r), sprintf('%s',ptxt)},'fontsize',8,'color','k','horizo','center','interpreter','latex')
           
           
           % least-squares fit
           pf = polyfit(exe,obs,1);
           x = linspace(-1,1,1000);
           y = polyval(pf,x);
           
           line([-1 1],[-1 1],'color',[1 1 1]*0.5,'linew',1.5,'linestyle','--'); % unity line
           line(x,y,'color',[1 1 1]*0,'linew',1.5)

           ax_inset.XTick = -1:1:1;
           ax_inset.YTick = -1:1:1;
           ax_inset.XAxis.TickLabelInterpreter='latex';
           ax_inset.YAxis.TickLabelInterpreter='latex';
           axis(ax_inset,[-1 1 -1 1]);
           axis square;
            
        end
            
        % shuffles
        for iter=1:niter
            [rr,~] = corr(exe,obs(randperm(length(obs))),'rows','complete','type',corrtype);
            Rshuff(iter,e,celltype) = rr;
        end
       
        R(celltype,e) = r;
        P(celltype,e) = p;
       
    end
end

hb  = plot(ax1,R');

% plot 5th and 95th percentiles of shuffled distribution
% hb5 = plot(ax1,squeeze(prctile(Rshuff,5,1)));
hb95 = plot(ax1,squeeze(prctile(Rshuff,95,1)));

clear leg
hold on
for i=1:length(hb)
    set(hb(i),'color',EpochAct.popncols(i,:),'linewidth',2,'marker','o');

%     set(hb5(i),'color',EpochAct.popncols(i,:),'linewidth',1.2,'linestyle',':');
    set(hb95(i),'color',EpochAct.popncols(i,:),'linewidth',1.2,'linestyle',':');
    
    text(ax1,6.5,1.05-(i-1)*0.1,sprintf('%s, n=%d',EpochAct.popnlabels{i},Nsel(i)),'color',EpochAct.popncols(i,:),'interpreter','latex')
end

ax1.XTick = 1:length(EpochAct.labs(epochs));
ax1.XTickLabel = EpochAct.labs(epochs);
ax1.XAxis.TickLabelInterpreter = 'latex';
ax1.YAxis.TickLabelInterpreter = 'latex';
ax1.YTick = -2:.2:2;
xtickangle(ax1,45);
axis(ax1,[0.75 numel(epochs)+0.25 -0.6 1.1])

ylabel(ax1,sprintf('%s Correlation Coefficient',corrtype),'interpreter','latex');
text(ax1,3,0.27,'Shuffles','interpreter','latex','fontsize',12)
tidyaxes(ax1);

% filename = 'Fig6a_CorrCoefLine';
% savefig(filename)

%% some statistics

for celltype=1:length(EpochAct.popnlabels)
    p = permtest(R(celltype,6),Rshuff(:,6,celltype),0.05,'upper',1);
%     p = permtest(R(celltype,6),max(Rshuff(:,6,:),[],3),0.05,'upper',1);
    fprintf('%s:\t p = %.3f\n',EpochAct.popnlabels{celltype},p)
end

%% scatter plots (Supplementary)

if plotsupps

clear R P C
epochs = [3 6 8]; % ObjCue, Early Reach, Hold
let = 'ABC';

figure('color','w','units','centimeters','position',[10 10 15 15]); clf
ha = tight_subplot(numel(epochs),length(EpochAct.popnlabels),[.1 .1],[.1 .1],[.1 .1]);

k=1;
for e=1:length(epochs)
    for celltype=1:length(EpochAct.popnlabels)
        
        selc = EpochAct.popn==celltype & ismirror;
        Nsel(celltype) = nnz(selc);
        
        exe1 = NormFRs(selc,epochs(e),CONDS.hm=='m'&CONDS.gng==1);
        obs1 = NormFRs(selc,epochs(e),CONDS.hm=='h'&CONDS.gng==1);   
        
        exe=exe1(:);
        obs=obs1(:);
        
        [r,p] = corr(exe,obs,'rows','complete','type',corrtype);

        axes(ha(k)); hold on
%         scatter(ha(k),exe,obs,12,'k','filled');
        scatter(ha(k),exe1(:,1),obs1(:,1),12,'r','filled'); % PG
        scatter(ha(k),exe1(:,2),obs1(:,2),12,'b','filled'); % WHG

        ptxt=sprintf('p = %.g',p);
        tt=text(1,-0.9,{sprintf('r = %.2f',r),ptxt},'color','k','horizo','right','interpreter','latex');
        if p > 0.05
            tt.Color=[1 1 1]*.5;
        end
        
        % least-squares fit
        pf = polyfit(exe,obs,1);
        x = linspace(-1,1,1000);
        y = polyval(pf,x);
        
        line([-1 1],[-1 1],'color',[1 1 1]*0.5,'linew',1.5,'linestyle','--'); % unity line
        line(x,y,'color',[1 1 1]*0,'linew',1.5)
        
        ha(k).XTick = -1:1:1;
        ha(k).YTick = -1:1:1;
        ha(k).XAxis.TickLabelInterpreter='latex';
        ha(k).YAxis.TickLabelInterpreter='latex';
        axis(ha(k),[-1 1 -1 1]);
        offsetAxes(ha(k));
        axis square;
        tidyaxes(ha(k));
        
        if celltype==1
            text(-1.25,1.15,let(ceil(k/length(EpochAct.popnlabels))),'fontsize',14,'fontweight','bold','verti','bottom')
        elseif celltype==3
            yl=ylabel(EpochAct.labs(epochs(e))); 
            yl.Interpreter = 'latex';
            yl.Position(1) = +1.5;
        end
        
        if e==1
            title(EpochAct.popnlabels{celltype},'interpreter','latex')
            if celltype==1
                text(-0.9,0.85,'PG','color','r','interpreter','latex','fontsize',10)
                text(-0.9,0.55,'WHG','color','b','interpreter','latex','fontsize',10)
            end
        end
        R(celltype,e) = r;
        P(celltype,e) = p;

        k=k+1;
    end
end
[ax,h]=suplabel('Net Normalised Execution Activity','x',[.1 .10 .84 .84]); set(h,'FontSize',14,'interpreter','latex')
[ax,h]=suplabel('Net Normalised Observation Activity','y',[.11 .10 .84 .84]); set(h,'FontSize',14,'interpreter','latex')

end

% filename = 'Fig6Supp_CorrCoefScatter';
% savefig(filename)