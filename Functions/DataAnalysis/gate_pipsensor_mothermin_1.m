function badtraces=gate_pipsensor(data,tracestats,motherstats,firstthreshold,maxthreshold,minthreshold,negthreshold)
%{
hist(max(signals,[],2),0:2:202); xlim([0 200]);
hist(max(signals,[],2),0:10:2010); xlim([0 2000]);
hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);

firstvals=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstvals(i)=signals(i,tracestats(i,1));
end
hist(firstvals,100);
%}

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(data,1);
signals=data;
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%totalmin=min(signals,[],2)<minthreshold & min(signals,[],2)>negthreshold;
%negmin=min(signals,[],2)>negthreshold;
totalmax=max(signals,[],2)>maxthreshold;
firstmax=ones(numtraces,1)*NaN;
mothermin=ones(numtraces,1)*NaN;
negmin=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstmax(i)=signals(i,tracestats(i,1))>firstthreshold;
    mothermin(i)=min(signals(i,motherstats(i,1):motherstats(i,2)))<minthreshold;
    negmin(i)=min(signals(i,motherstats(i,1):tracestats(i,2)))>negthreshold;
end
goodtraces= negmin & totalmax & firstmax & mothermin;
badtraces=~goodtraces;