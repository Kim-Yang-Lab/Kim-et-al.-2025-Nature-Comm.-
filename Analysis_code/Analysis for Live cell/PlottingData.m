%% load Data
clear;
root=''; % Path
resultdir=[root 'Results/'];
load([resultdir 'Data_' num2str(option) '.mat']);
numCondi=length(sensor);
%% checking data
condi=10;
figure;
subplot(1,2,1);
plot(sensor(condi).APC_smooth','b');
subplot(1,2,2);
plot(sensor(condi).Cdt1_smooth','b');
%% checking data
condi=1;
figure; histogram(sensor(condi).APC_smooth(:,1));
%% Selecting G2-phase cells based on APC and Cdt1
numCondi=length(sensor);
for condi=1:numCondi;
    cellNum=size(sensor(condi).APC_smooth,1);
    for i=1:cellNum
        idxLast=find(sum(~isnan(sensor(condi).APC_smooth(i,:)),1) > 0, 1 ,'last');
        APC{condi}(i)=nanmean(sensor(condi).APC_smooth(i,45:50));
        APC_end{condi}(i)=sensor(condi).APC_smooth(i,idxLast);
        Cdt1_before{condi}(i)=nanmean(sensor(condi).norm_Cdt1_smooth(i,1:10));
        Cdt1_after{condi}(i)=nanmean(sensor(condi).norm_Cdt1_smooth(i,35:45));
    end
end
%% Selecting G2-phase cells based on APC and Cdt1
conti=1;
thresh_APC=100;
thresh_Cdt1=0.4;
%nucArea=nanmean(sensor(condi).nucArea(:,1:20),2);
figure; subplot(1,3,1);
histogram(APC{condi},100); vline(thresh_APC);
subplot(1,3,2);
histogram(Cdt1_before{condi}); vline(thresh_Cdt1); xlim([0 1]);
subplot(1,3,3);
histogram(Cdt1_after{condi}); vline(thresh_Cdt1); xlim([0 1]);
%% Selecting G2-phase cells based on APC and Cdt1
thresh_APC=100;
figure; [h]=histogram(APC{condi},100,'Normalization','probability','BinLimits',[0,1000]);
%xlim([80,700])
vline([thresh_APC ]);
ylabel('Cells'); xlabel('APC signal');
lowNucArea=sum(h.Values(h.BinEdges>=h.BinLimits(1) & h.BinEdges<thresh_APC));
highNucArea=sum(h.Values(h.BinEdges>=thresh_APC & h.BinEdges<h.BinLimits(2)));
%title(['Low area: ' num2str(lowNucArea*100) '% Middle area: ' num2str(middleNucArea*100) '%' '% High area: ' num2str(highNucArea*100) '%'])
%%
for condi=1:numCondi;
    gating{condi}=APC{condi}>thresh_APC & Cdt1_before{condi}<thresh_Cdt1 & Cdt1_before{condi}>0 & Cdt1_after{condi}>0.1;
end
%%
condi=5;
numberMitosis=sensor(condi).numMitosis;
gate_0=sum(sensor(condi).norm_Cdt1_smooth(:,10:40)>0.35,2)>1 | sum(sensor(condi).norm_Cdt1_smooth(:,80:90)<0.2,2)>1;
gate_1=~gate_0 & gating{condi}' & nanmean(sensor(condi).APC_smooth(:,38:42),2)>600 & (numberMitosis>=0);
figure;
subplot(1,2,1);
plot(sensor(condi).APC_smooth(gate_1,:)','k');
xlim([0 190]); ylim([0 10000]); vline(40); vline(160);
ylabel('APC/C-degron intensity'); xlabel('Timse since mitogen stimulation');
set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
subplot(1,2,2); hold on;
plot(sensor(condi).norm_Cdt1_smooth(gate_1,:)','k');
xlim([0 190]); ylim([0 1]); vline(40); vline(160);
ylabel('Cdt1-degron intensity'); xlabel('Timse since mitogen stimulation');
set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
%% Plot single-cell traces classified based on APC/C reactivation (Threshold needs to be adjusted based on the control condition)
clear gate_1 gate_2
close all;
count=0;
for condi=[1:5];
    count=count+1;
    clear Idx; %close all;
    cellNum=50;
    color=jet(cellNum);
    numberMitosis=sensor(condi).numMitosis;
    gate_0=sum(sensor(condi).norm_Cdt1_smooth(:,5:15)>0.35,2)>1 | sum(sensor(condi).norm_Cdt1_smooth(:,80:90)<0.2,2)>1;
    gate_0_1=sensor(condi).mitosisFrame(:,1)<45;
    if option==1
        gate_1{condi}=~gate_0 & ~gate_0_1 & gating{condi}' & nanmean(sensor(condi).APC_smooth(:,38:42),2)>1500 & (numberMitosis>=0);
        numCell_gate1=sum(gate_1{condi});
        numTotal=numCell_gate1;
        Idx_gate1=find(gate_1{condi});
        MitoFrame{condi}=sensor(condi).mitosisFrame(Idx_gate1,1)-40;
        plotCells_gate1=randsample(Idx_gate1,cellNum);
        figure; subplot(1,2,1); hold on
        for j=1:cellNum
            plot(sensor(condi).APC_smooth(plotCells_gate1(j),:)','k');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate1(j),:)));
            if numMito
                for k=1%:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate1(j),k),sensor(condi).APC_smooth(plotCells_gate1(j),sensor(condi).mitosisFrame(plotCells_gate1(j),k)),'r','filled');
                end
            end
        end
        xlim([0 190]); ylim([0 6000]); vline([38 40]);
        ylabel('APC/C-degron intensity'); xlabel('Timse (hr)');
        set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
        %title(['Red Cells: ' num2str(ratio_gate1*100)]);
        subplot(1,2,2); hold on
        for j=1:cellNum
            plot(sensor(condi).norm_Cdt1_smooth(plotCells_gate1(j),:)','k');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate1(j),:)));
            if numMito
                for k=1%:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate1(j),k),sensor(condi).norm_Cdt1_smooth(plotCells_gate1(j),sensor(condi).mitosisFrame(plotCells_gate1(j),k)),'r','filled');
                end
            end
        end
        ylabel('Cdt1-degron intensity'); xlabel('Timse since mitogen stimulation');
        set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
        xlim([0 190]); ylim([0 1]); vline([38 40]);
        figu(0.3,0.9)
        %title(['Black Cells: ' num2str(ratio_gate2*100)]);
        title(conditions{condi});
        figu(0.3,0.9)
    elseif option==2
        gate_1{condi}=~gate_0 & ~gate_0_1 & gating{condi}' & nanmean(sensor(condi).APC_smooth(:,38:42),2)>1500 & nanmean(sensor(condi).APC_smooth(:,158:160),2)<3500 & (numberMitosis>=0);
        gate_2{condi}=~gate_0 & ~gate_0_1 & gating{condi}' & nanmean(sensor(condi).APC_smooth(:,38:42),2)>1500 & nanmean(sensor(condi).APC_smooth(:,118:122),2)>3000 & nanmean(sensor(condi).APC_smooth(:,158:160),2)>5000 & sensor(condi).norm_Cdt1_smooth(:,120)>0.35 & (numberMitosis>=0);
        numCell_gate1=sum(gate_1{condi});
        numCell_gate2=sum(gate_2{condi});
        numTotal=numCell_gate1+numCell_gate2;
        ratio_gate1=numCell_gate1/numTotal;
        ratio_gate2=numCell_gate2/numTotal;
        num_gate1=round(ratio_gate1*cellNum);
        num_gate2=round(ratio_gate2*cellNum);
        %num_gate1=cellNum;
        %num_gate2=cellNum;
        Idx_gate1=find(gate_1{condi});
        Idx_gate2=find(gate_2{condi});
        plotCells_gate1=randsample(Idx_gate1,num_gate1);
        plotCells_gate2=randsample(Idx_gate2,num_gate2);
        APCreactivation(count)=ratio_gate1*100;
        figure; subplot(1,2,1); hold on
        for j=1:num_gate1
            plot(sensor(condi).APC_smooth(plotCells_gate1(j),:)','r');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate1(j),:)));
            if numMito
                for k=1:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate1(j),k),sensor(condi).APC_smooth(plotCells_gate1(j),sensor(condi).mitosisFrame(plotCells_gate1(j),k)),'r','filled');
                end
            end
        end
        xlim([0 160]); ylim([0 10000]); vline(40); vline(160);
        ylabel('APC/C-degron intensity'); xlabel('Timse (hr)');
        set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
        title(['Red Cells: ' num2str(ratio_gate1*100)]);
        for j=1:num_gate2
            plot(sensor(condi).APC_smooth(plotCells_gate2(j),:)','k');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate2(j),:)));
            if numMito
                for k=1%:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate2(j),k),sensor(condi).APC_smooth(plotCells_gate2(j),sensor(condi).mitosisFrame(plotCells_gate2(j),k)),'k','filled');
                end
            end
        end
        subplot(1,2,2); hold on
        for j=1:num_gate1
            plot(sensor(condi).norm_Cdt1_smooth(plotCells_gate1(j),:)','r');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate1(j),:)));
            if numMito
                for k=1:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate1(j),k),sensor(condi).norm_Cdt1_smooth(plotCells_gate1(j),sensor(condi).mitosisFrame(plotCells_gate1(j),k)),'r','filled');
                end
            end
        end
        for j=1:num_gate2
            plot(sensor(condi).norm_Cdt1_smooth(plotCells_gate2(j),:)','k');
            numMito=sum(~isnan(sensor(condi).mitosisFrame(plotCells_gate2(j),:)));
            if numMito
                for k=1%:numMito
                    scatter(sensor(condi).mitosisFrame(plotCells_gate2(j),k),sensor(condi).norm_Cdt1_smooth(plotCells_gate2(j),sensor(condi).mitosisFrame(plotCells_gate2(j),k)),'k','filled');
                end
            end
        end
        ylabel('Cdt1-degron intensity'); xlabel('Timse since mitogen stimulation');
        set(gca,'XTick',[0:40:300]); xPos=[0:8:100]; set(gca,'XTickLabel',xPos);
        xlim([0 160]); ylim([0 1]); vline(40); vline(160);
        figu(0.3,0.9)
        title(['Black Cells: ' num2str(ratio_gate2*100)]);
        title(conditions{condi});
        figu(0.3,0.9)
    end
end
%save([resultdir 'results_20230830.mat'],'mitosisTime','mitosisPrc','APCreactivation');
