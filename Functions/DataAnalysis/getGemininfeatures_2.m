function [risetime,badtraces]=getGemininfeatures_1(sampletraces,samplestats)
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=sigstore;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
    %signal=signal(firstframe(i):lastframe(i));
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prebuffer=5;
    badtraces(i)=signal(prebuffer+1)>0.02;
    
    %%%%%% general requirements
    postbuffer=10;
    gate=zeros(1,numframes);
    for j=prebuffer:numframes-postbuffer
        presentheight=signal(j)<0.2;
        minfutureheight=min(signal(j+1:end))>signal(j);
        bufferfutureheight=signal(j+postbuffer)>signal(j)+0.5;
        maxfutureheight=max(signal(j+1:end))>(signal(j)+1);
        gate(j)=presentheight & minfutureheight & bufferfutureheight & maxfutureheight;
    end
    
    %sig_revslope=10*getslope_reverse(signal,1:10);
    %sig_fwdslope=10*getslope_forward_avg(signal,1:5);
    sig_fwdslope=10*getslope_forward_avg(signal,1:postbuffer);
    sig_time=(1:length(signal))/tracelength;
    filterscore=2*sig_fwdslope-4*signal+sig_time+6;
    filterscore=filterscore.*gate;
    altstore(i,firstframe(i):lastframe(i))=sig_fwdslope;
    %tempsearch=find(sig_fwdslope>0.05 & abs(signal)<0.03,1,'last');
    %tempsearch=find(filterscore>0.05,1,'first');
    filtermax=max(filterscore);
    tempsearch=find(filterscore==filtermax,1,'first');
    %if isempty(tempsearch) || signal(end)<0.05 %0.05
    if isempty(tempsearch)
        risetime(i)=NaN;
    else
        risetimedb(i)=tempsearch;
        risetime(i)=tempsearch+firstframe(i)-1; %return absolute POI rather than relative to mitosis
    end
end
%keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    ylim([0 1]);
    if badtraces(i)==1 || isnan(risetime(i))
        continue;
    end
    hold on;
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
    %plot(1:length(frames),altstore(i,frames),'r');
end
%}