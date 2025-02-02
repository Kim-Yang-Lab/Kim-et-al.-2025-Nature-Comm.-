function [minval,risetime,riseslope,badtraces]=getCdk2features(sampletraces)
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
minval=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
riseslope=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=find(~isnan(signal),1,'first');
    lastframe(i)=find(~isnan(signal),1,'last');
    signal=smooth(signal(firstframe(i):lastframe(i)))';
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% calc risetime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig_fwdslope=getslope_forward_avg(signal,6:10);
    %sig_fwdslope_long=getslope_forward_avg(signal,11:20);
    sig_height=signal;
    sig_time=(1:length(signal))/length(signal);
    heightgate=sig_height<(min(sig_height)+0.3*range(sig_height));
    filter=2*sig_fwdslope-sig_height+sig_time+2;
    %filter=sig_fwdslope+sig_fwdslope_long-2*sig_height+sig_time+2;
    filter=filter.*heightgate;
    %filter=4*sig_fwdslope+2*sig_fwdslope_long-sig_height+100;
    filter(end-9:end)=min(filter);
    filter(1)=min(filter); %remove noisy first signal
    risetime(i)=find(filter==max(filter),1,'first');
    %%% calc minval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %minval(i)=signal(risetime(i));
    minval(i)=min(signal(1:risetime(i)));
    if minval(i)>1
        badtraces(i)=1;
        continue;
    end
    %%% calc riseslope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    riseslope(i)=sig_fwdslope(risetime(i));
    %filter=filter-min(filter);
    %filter=filter/max(filter);
    altstore1(i,firstframe(i):lastframe(i))=filter;
end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(frames,sigstore(i,frames));
    axis([frames(1) frames(end) 0 2]);
    if badtraces(i)==1
        continue;
    end
    hold on;
    %plot(frames,altstore1(i,frames),'r');
    %plot(risetime(i),minval(i),'go','markerfacecolor','g','markersize',6);
    plot(risetime(i):risetime(i)+9,minval(i)+riseslope(i)*(0:9),'r');
end


plot(x,sigstore(22,:),'linewidth',4); xlim([1 length(signal)]);;
hold on;
plot(degstarts(22),sigstore(22,degstarts(22)),'go','markerfacecolor','g','markersize',12);
plot(degends(22),sigstore(22,degends(22)),'ro','markerfacecolor','r','markersize',12);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');

good=find(badtraces==0 & degends<length(signal));
bad=find(badtraces==0 & degends==length(signal));
goodvals=ones(numel(good),1)*NaN;
for i=1:numel(good)
    g=good(i);
    goodvals(i)=altstore(g,degends(g));
end
figure,hist(goodvals);
badvals=ones(numel(bad),1)*NaN;
for i=1:numel(bad)
    b=bad(i);
    badvals(i)=altstore(b,degends(b));
end
figure,hist(badvals);
%}