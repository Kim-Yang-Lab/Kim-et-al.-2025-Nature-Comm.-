function comparetrends_1(datatotal,xvals,setb,setr,xstring,ystring,ylimits,smoothoption)
datab=datatotal(setb,:);
datar=datatotal(setr,:);
figure, hold on;
for i=1:size(datab,1)
    plotdatab=datab(i,:);
    if smoothoption
        realframes=find(~isnan(datab(i,:)));
        plotdatab(realframes)=smooth(datab(i,realframes));
    end
    line(xvals,plotdatab,'color','b');
end
for i=1:size(datar,1)
    plotdatar=datar(i,:);
    if smoothoption
        realframes=find(~isnan(datar(i,:)));
        plotdatar(realframes)=smooth(datar(i,realframes));
    end
    line(xvals,plotdatar,'color','r');
end
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
%saveas(gcf,'h:\Downloads\FigCompareLines.jpg');

figure, hold on;
binstep=0.2; %absolute bin width (usu 1frame=0.2hrs)
bincurveshade(xvals,datab,binstep,'b');
bincurveshade(xvals,datar,binstep,'r');
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
%saveas(gcf,'h:\Downloads\FigTrends.jpg');