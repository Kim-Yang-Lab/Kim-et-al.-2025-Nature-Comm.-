function Data=CellCycleGating(conditions,datadir,option,Celloption)
showylabel=1;
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
bmin=5; bmax=13; %pRb
%bmin=0; bmax=1; %CDK2
bstep=(bmax-bmin)/35; %pRb:50 pRb/tRb:100
bin=bmin:bstep:bmax;
[~,idx1]=min(abs(bin-1));
numbins=numel(bin);
binfill=[bin fliplr(bin)];
namecount=cell(uniquecondnum,1);
colorcode='krgcmb';
colors=colorcode;
fontsizevar=10; %boxplots:12 histogram:12
if showylabel
    mlvar=0.25; %default 0.25
    %pphvar=2; %4col per slide
    %pphvar=2.5;
    pphvar=4; %2col per slide
    %phvar=6; %1col per slide
else
    mlvar=0.1;
    %pphvar=1.75; %4col per slide
    %pphvar=2;
    pphvar=3.5; %2col per slide
end
pmat=cell(uniquecondnum,1);
boxplotdata=[];
offset=[];
for i=1:uniquecondnum
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    for i=condrow'
        rowmat=cell2mat(conditions(i,2));
        colmat=cell2mat(conditions(i,3));
        sitemat=cell2mat(conditions(i,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];

                    %[rowstring,colstring,sitestring]=wellnum2strRCS_3(row,col,site);
                    %shot=[rowstring,colstring,'_',sitestring];
                    %load([datadir,shot,'.mat'],'IFdata');
                    if exist([datadir,'IF_',shot,'.mat'])
                        %load([datadir,shot,'.mat'],'IFdata');
                        load([datadir,'IF_',shot,'.mat'],'IFdata');
                        Alldata=[Alldata;IFdata];
                    else
                        continue
                    end
                end
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    %Hoechstval=Alldata(:,3);
    %histogram(Hoechstval)
    %EdUval=Alldata(:,4); APC
    %histogram(EdUval);
    %%%%%%%%%%%%%%%%%%change this%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    EdUval=Alldata(:,6);

    EdUval(EdUval<1)=1;
    EdUval=log2(EdUval);
    %histogram(EdUval);
    %dscatter(Hoechstval,EdUval)
    %%% Varying cell cycle gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plotting boxes based on picture (change this)
   
        if option==1
            G1minH=[000000,000000,000000,000000,000000,000000];
            G1maxH=[5500000,5500000,5500000,5500000,5500000,5500000];
            G1minE=[0,0,0,0,0,0]; G1maxE=[10,10,10,10,10,10];
            SminH=[000000,000000,000000,000000,000000,000000];
            SmaxH=[10000000,10000000,10000000,10000000,10000000,10000000];
            SminE=[10,10,10,10,10,10];
            SmaxE=[16,16,16,16,16,16];
            G2minH=[5500000,5500000,5500000,5500000,5500000,5500000];
            G2maxH=[10000000,10000000,10000000,10000000,10000000,10000000];
            G2minE=[0,0,0,0,0,0]; G2maxE=[10,10,10,10,10,10];
            xylim=[000000 10000000 0 16];
        elseif option==2
            G1minH=[000000,000000,000000,000000,000000,000000];
            G1maxH=[3500000,4000000,3500000,3500000,3600000,3500000];
            G1minE=[0,0,0,0,0,0]; G1maxE=[9,9,9,9,9,9];
            SminH=[000000,000000,000000,000000,000000,000000];
            SmaxH=[6400000,6400000,6400000,6400000,6400000,6400000];
            SminE=[9,9,9,9,9,9];
            SmaxE=[15,15,15,15,15,15];
            G2minH=[3500000,4000000,3500000,3500000,3600000,3500000];
            G2maxH=[6400000,6400000,6400000,6400000,6400000,6400000];
            G2minE=[0,0,0,0,0,0]; G2maxE=[9,9,9,9,9,9];
            xylim=[000000 6400000 0 15];
        end     

    %%% Report Cell-cycle percentages
    G1cells=Hoechstval>G1minH(i) & Hoechstval<G1maxH(i) & EdUval>G1minE(i) & EdUval<G1maxE(i);
    Scells=Hoechstval>SminH(i) & Hoechstval<SmaxH(i) & EdUval>SminE(i) & EdUval<SmaxE(i);
    G2cells=Hoechstval>G2minH(i) & Hoechstval<G2maxH(i) & EdUval>G2minE(i) & EdUval<G2maxE(i);

    numG1=sum(G1cells);
    numS=sum(Scells);
    numG2=sum(G2cells);
    numall=numG1+numS+numG2;
    percG1=100*numG1/numall;
    percS=100*numS/numall;
    percG2=100*numG2/numall;
    fprintf('G1 = %0.0f\n',percG1);
    fprintf('S = %0.0f\n',percS);
    fprintf('G2 = %0.0f\n',percG2);

    disp(['G1 = ' num2str(percG1)]);
    Data{i}(1)=percG1;
    Data{i}(2)=percS;
    Data{i}(3)=percG2;

    %     if ismember(option,[1 2 3 5 6])
    %         Hoechst_plot=Hoechstval;
    %         EdU_plot=EdUval;
    %         subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.2,'SH',0.03); %horizontal
    %         %subaxis(uniquecondnum,1,i,'ML',0.3,'MR',0.03,'MT',0.1,'MB',0.2,'SV',0.1); %vertical
    %         dscatter_gray(Hoechst_plot,EdU_plot);
    %     else
    %         if ismember(option,[8])
    %             numCell=1000;
    %         else
    numCell=1500;
    %         end
    idx_G1=find(G1cells); gating_G1=randsample(idx_G1,round(numCell*percG1/100));
    idx_S=find(Scells); gating_S=randsample(idx_S,round(numCell*percS/100));
    idx_G2=find(G2cells); gating_G2=randsample(idx_G2,round(numCell*percG2/100));
    gating=[gating_G1; gating_S; gating_G2];
    Hoechst_plot=Hoechstval(gating);
    EdU_plot=EdUval(gating);
    subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.2,'SH',0.03); %horizontal
    %subaxis(uniquecondnum,1,i,'ML',0.3,'MR',0.03,'MT',0.1,'MB',0.2,'SV',0.1); %vertical
    dscatter_gray(Hoechst_plot,EdU_plot);
    %     end
    %     else
    %         subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.2,'SH',0.03); %horizontal
    %         %subaxis(uniquecondnum,1,i,'ML',0.3,'MR',0.03,'MT',0.1,'MB',0.2,'SV',0.1); %vertical
    %         dscatter_gray(Hoechstval,EdUval);
    %     end

    %     Hoechst_plot=Hoechstval(gating);
    %     EdU_plot=EdUval(gating);
    %     subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.2,'SH',0.03); %horizontal
    %     %subaxis(uniquecondnum,1,i,'ML',0.3,'MR',0.03,'MT',0.1,'MB',0.2,'SV',0.1); %vertical
    %     dscatter_gray(Hoechst_plot,EdU_plot);
    %hist(Hoechstval,-50:50:1550); xlim([0 1500]);
    %%% draw cell cycle gates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on;
    rectangle('Position',[G1minH(i),G1minE(i),G1maxH(i)-G1minH(i),G1maxE(i)-G1minE(i)],'EdgeColor','r','linewidth',1);
    rectangle('Position',[SminH(i),SminE(i),SmaxH(i)-SminH(i),SmaxE(i)-SminE(i)],'EdgeColor','r','linewidth',1);
    rectangle('Position',[G2minH(i),G2minE(i),G2maxH(i)-G2minH(i),G2maxE(i)-G2minE(i)],'EdgeColor','r','linewidth',1);
    axis(xylim);

    %%% horizontal
    %title(char(uniquenames(i)),'fontsize',fontsizevar);
    %%% vertical
    y=ylabel(char(uniquenames(i)),'rot',0,'fontsize',fontsizevar);
    set(y,'Units','Normalized','Position',[0.55,1,0]);
    %%%
    set(gca,'fontsize',fontsizevar);
    namecount{i}=char(uniquenames(i));
end

%set(gcf,'color','w','PaperPosition',[0 0 8 4]); %horizontal
set(gcf,'color','w','PaperPosition',[0 0 6 8]); %vertical

%textpos=[0.15 0.55 0.3 0.15];
%annotation('textbox',textpos,'String',['p-value = ',num2str(p)]);
% saveas(gcf,'D:\Downloads\Fig.jpg');
%print('-depsc','-tiff','-r300','h:\Downloads\Fig.eps');
%%%
end